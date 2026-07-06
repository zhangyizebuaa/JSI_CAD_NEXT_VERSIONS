#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GeomAPI_IntSS.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <OSD_ThreadPool.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "nurbs_data.hpp"

namespace {

struct KnotData {
  std::vector<double> values;
  std::vector<int> multiplicities;
};

KnotData compress_knots(const float *knots, int count) {
  KnotData result;
  for (int i = 0; i < count; ++i) {
    const double value = knots[i];
    if (result.values.empty() || value != result.values.back()) {
      result.values.push_back(value);
      result.multiplicities.push_back(1);
    } else {
      ++result.multiplicities.back();
    }
  }
  return result;
}

template <typename T>
NCollection_Array1<T> make_array1(const std::vector<T> &values) {
  NCollection_Array1<T> result(1, static_cast<int>(values.size()));
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    result.SetValue(i + 1, values[static_cast<size_t>(i)]);
  }
  return result;
}

occ::handle<Geom_BSplineSurface> make_surface(const NurbsFace &face) {
  const int flags = face.ivec[0];
  const int u_degree = face.ivec[1];
  const int u_poles = face.ivec[2];
  const int u_knot_count = face.ivec[3];
  const int v_degree = face.ivec[4];
  const int v_poles = face.ivec[5];
  const int v_knot_count = face.ivec[6];
  const size_t xyz_offset = static_cast<size_t>(u_knot_count + v_knot_count);
  const size_t weight_offset = xyz_offset + static_cast<size_t>(u_poles) * v_poles * 3;

  if (u_degree <= 0 || v_degree <= 0 || u_poles <= 1 || v_poles <= 1 ||
      weight_offset > static_cast<size_t>(face.ndata)) {
    throw std::runtime_error("invalid JSI NURBS data");
  }

  const KnotData u_knots = compress_knots(face.data, u_knot_count);
  const KnotData v_knots = compress_knots(face.data + u_knot_count, v_knot_count);
  const auto occt_u_knots = make_array1(u_knots.values);
  const auto occt_v_knots = make_array1(v_knots.values);
  const auto occt_u_mults = make_array1(u_knots.multiplicities);
  const auto occt_v_mults = make_array1(v_knots.multiplicities);

  NCollection_Array2<gp_Pnt> poles(1, u_poles, 1, v_poles);
  const float *control_points = face.data + xyz_offset;
  for (int v = 0; v < v_poles; ++v) {
    for (int u = 0; u < u_poles; ++u) {
      const size_t index = static_cast<size_t>(u + u_poles * v) * 3;
      poles.SetValue(u + 1, v + 1,
                     gp_Pnt(control_points[index], control_points[index + 1], control_points[index + 2]));
    }
  }

  if ((flags & 2) == 0) {
    return new Geom_BSplineSurface(poles, occt_u_knots, occt_v_knots, occt_u_mults, occt_v_mults, u_degree,
                                   v_degree);
  }

  if (weight_offset + static_cast<size_t>(u_poles) * v_poles > static_cast<size_t>(face.ndata)) {
    throw std::runtime_error("invalid rational JSI NURBS data");
  }
  NCollection_Array2<double> weights(1, u_poles, 1, v_poles);
  const float *raw_weights = face.data + weight_offset;
  for (int v = 0; v < v_poles; ++v) {
    for (int u = 0; u < u_poles; ++u) {
      weights.SetValue(u + 1, v + 1, raw_weights[u + u_poles * v]);
    }
  }
  return new Geom_BSplineSurface(poles, weights, occt_u_knots, occt_v_knots, occt_u_mults, occt_v_mults,
                                 u_degree, v_degree);
}

double mean(const std::vector<double> &values) {
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const size_t middle = values.size() / 2;
  return values.size() % 2 ? values[middle] : (values[middle - 1] + values[middle]) * 0.5;
}

template <typename F>
std::vector<double> measure_ms(int repeats, F &&operation) {
  std::vector<double> samples;
  samples.reserve(static_cast<size_t>(repeats));
  for (int i = 0; i < repeats; ++i) {
    const auto begin = std::chrono::steady_clock::now();
    operation();
    const auto end = std::chrono::steady_clock::now();
    samples.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
  }
  return samples;
}

std::pair<NurbsFace *, NurbsFace *> distance_case(int case_id) {
  NurbsFace *side1[] = {&dist_case1_obj1, &dist_case2_obj1, &dist_case_hard_parallel_obj1,
                        &dist_hard_wave_obj1};
  NurbsFace *side2[] = {&dist_case1_obj2, &dist_case2_obj2, &dist_case_hard_parallel_obj2,
                        &dist_hard_wave_obj2};
  if (case_id < 1 || case_id > 4) throw std::runtime_error("distance case must be 1..4");
  return {side1[case_id - 1], side2[case_id - 1]};
}

std::pair<NurbsFace *, NurbsFace *> intersection_case(int case_id) {
  NurbsFace *side1[] = {&inter_case1_obj1, &inter_case2_obj1};
  NurbsFace *side2[] = {&inter_case1_obj2, &inter_case2_obj2};
  if (case_id < 1 || case_id > 2) throw std::runtime_error("intersection case must be 1..2");
  return {side1[case_id - 1], side2[case_id - 1]};
}

void print_samples(const std::vector<double> &samples) {
  std::cout << "samples_ms=";
  for (size_t i = 0; i < samples.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << samples[i];
  }
  std::cout << "\nmean_ms=" << mean(samples) << "\nmedian_ms=" << median(samples) << '\n';
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
              << " distance|intersection|section CASE [REPEATS] [THREADS]\n";
    return 2;
  }

  try {
    const std::string mode = argv[1];
    const int case_id = std::stoi(argv[2]);
    const int repeats = argc >= 4 ? std::max(1, std::stoi(argv[3])) : 5;
    const int threads = argc == 5 ? std::max(1, std::stoi(argv[4])) : 1;
    OSD_ThreadPool::DefaultPool(threads)->SetNbDefaultThreadsToLaunch(threads);
    const auto faces = mode == "distance" ? distance_case(case_id) : intersection_case(case_id);
    const bool parallel = threads > 1 && mode != "intersection";

    occ::handle<Geom_BSplineSurface> surface1;
    occ::handle<Geom_BSplineSurface> surface2;
    TopoDS_Face face1;
    TopoDS_Face face2;
    const auto construction = measure_ms(1, [&] {
      surface1 = make_surface(*faces.first);
      surface2 = make_surface(*faces.second);
      face1 = BRepBuilderAPI_MakeFace(surface1, faces.first->get_ubegin(), faces.first->get_uend(),
                                      faces.first->get_vbegin(), faces.first->get_vend(), 1.0e-7);
      face2 = BRepBuilderAPI_MakeFace(surface2, faces.second->get_ubegin(), faces.second->get_uend(),
                                      faces.second->get_vbegin(), faces.second->get_vend(), 1.0e-7);
    });

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "engine=OpenCASCADE\nmode=" << mode << "\ncase=" << case_id << "\nrepeats=" << repeats
              << "\nthreads=" << threads << "\nparallel=" << (parallel ? "true" : "false")
              << "\nconstruction_ms=" << construction.front() << '\n';

    if (mode == "distance") {
      double distance = 0.0;
      auto run = [&] {
        BRepExtrema_DistShapeShape extrema;
        extrema.LoadS1(face1);
        extrema.LoadS2(face2);
        extrema.SetMultiThread(threads > 1);
        if (!extrema.Perform()) throw std::runtime_error("OpenCASCADE distance failed");
        distance = extrema.Value();
      };
      run();
      const auto samples = measure_ms(repeats, run);
      std::cout << "distance=" << distance << '\n';
      print_samples(samples);
    } else if (mode == "intersection") {
      int curve_count = 0;
      auto run = [&] {
        GeomAPI_IntSS intersection(surface1, surface2, 1.0e-7);
        if (!intersection.IsDone()) throw std::runtime_error("OpenCASCADE intersection failed");
        curve_count = intersection.NbLines();
      };
      run();
      const auto samples = measure_ms(repeats, run);
      std::cout << "intersection_curves=" << curve_count << '\n';
      print_samples(samples);
    } else if (mode == "section") {
      int curve_count = 0;
      auto run = [&] {
        BRepAlgoAPI_Section intersection(face1, face2, false);
        intersection.SetRunParallel(threads > 1);
        intersection.Build();
        if (!intersection.IsDone()) throw std::runtime_error("OpenCASCADE intersection failed");
        curve_count = 0;
        for (TopExp_Explorer edge(intersection.Shape(), TopAbs_EDGE); edge.More(); edge.Next()) {
          ++curve_count;
        }
      };
      run();
      const auto samples = measure_ms(repeats, run);
      std::cout << "intersection_curves=" << curve_count << '\n';
      print_samples(samples);
    } else {
      throw std::runtime_error("mode must be distance, intersection, or section");
    }
  } catch (const Standard_Failure &failure) {
    std::cerr << "OpenCASCADE error: " << failure.what() << '\n';
    return 1;
  } catch (const std::exception &error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
