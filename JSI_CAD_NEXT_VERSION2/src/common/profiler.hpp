#ifndef CAD_PROFILER_HPP
#define CAD_PROFILER_HPP

#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace cad {
namespace profile {

inline bool Enabled() {
  static int enabled = []() -> int {
    const char* env = std::getenv("CAD_PROFILE");
    return (env != nullptr && env[0] != '\0' && env[0] != '0') ? 1 : 0;
  }();
  return enabled == 1;
}

class ScopeTimer {
 public:
  explicit ScopeTimer(const char* label)
      : label_(label), begin_(std::chrono::steady_clock::now()), enabled_(Enabled()) {}

  ~ScopeTimer() {
    if (!enabled_) return;
    const auto end = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(end - begin_).count();
    std::fprintf(stderr, "[cad-profile] %s: %.6f ms\n", label_, ms);
  }

 private:
  const char* label_;
  std::chrono::steady_clock::time_point begin_;
  bool enabled_;
};

}  // namespace profile
}  // namespace cad

#define CAD_PROFILE_SCOPE(name_literal) \
  ::cad::profile::ScopeTimer cad_profile_scope_timer_##__LINE__(name_literal)

#endif
