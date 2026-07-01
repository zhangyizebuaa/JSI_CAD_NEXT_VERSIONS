#include <iostream>
#include <chrono>

#include "nurbs_data.hpp"
#include "egads.h"
#include "egadsClasses.h"

NurbsFace nurbs_convert(const std::string& filename) {
    ego context, model;
    int status, nface, nbody;
    int oclass, mtype, *senses;
    ego geom, *faces, *bodies;
    NurbsFace resultFace;

    static char *classType[27] = {"CONTEXT", "TRANSFORM", "TESSELLATION",
                                "NIL", "EMPTY", "REFERENCE", "", "",
                                "", "", "PCURVE", "CURVE", "SURFACE", "", 
                                "", "", "", "", "", "", "NODE",
                                "EGDE", "LOOP", "FACE", "SHELL",
                                "BODY", "MODEL"};
    static char *surfType[11] = {"Plane", "Spherical", "Cylinder", "Revolution",
                               "Toroidal", "Trimmed" , "Bezier", "BSpline", 
                               "Offset", "Conical", "Extrusion"};

    status = EG_open(&context);
    if (status != EGADS_SUCCESS) {
        printf("EGADS open failed: %d\n", status);
        throw std::runtime_error("Failed to open EGADS context.");
    }

    status = EG_loadModel(context, 0, filename.c_str(), &model);
    if (status != EGADS_SUCCESS) {
        printf("EGADS load failure: %d\n", status);
        EG_deleteObject(model);
        EG_close(context);
        throw std::runtime_error("Failed to load model.");
    }
    
    status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies, &senses);
    if (status != EGADS_SUCCESS) {
        printf("EG_getTopology failed: %d\n", status);
        throw std::runtime_error("Failed to get model topology.");
    }

    for (int i = 0; i < nbody; i++) {
        status = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
        if (status != EGADS_SUCCESS || nface <= 0) {
            continue;
        }

        ego geom1, *bodies1;
        int nbody1;
        int oclass1, mtype1, *senses1;
        status = EG_getTopology(faces[0], &geom1, &oclass1, &mtype1, NULL, &nbody1, &bodies1, &senses1);
        if (status != EGADS_SUCCESS) {
            printf("EG_getTopology failed: %d\n", status);
            throw std::runtime_error("Failed to get model topology.");
        }
        ego geom2;
        int periodic, oclass2, mtype2, *ivec2;
        double limits[4], *rvec2;
        status = EG_getGeometry(geom1, &oclass2, &mtype2, &geom2, &ivec2, &rvec2);
        if (oclass2 == SURFACE && mtype2 == BSPLINE) {
            egadsSurface *lgeom = (egadsSurface *) geom1->blind;
            int len = lgeom->dataLen;
            int *ivec = lgeom->header;
            double *gdata = lgeom->data;

            // Populate NurbsFace struct
            for (int j = 0; j < 7; j++) {
                resultFace.ivec[j] = ivec[j];
            }
            resultFace.ndata = len;
            resultFace.data = new float[len];
            for (int j = 0; j < len; j++) {
                resultFace.data[j] = static_cast<float>(gdata[j]);
            }
            break;
        }
    }

    EG_deleteObject(model);
    EG_close(context);

    return resultFace;
}

int main(int argc, char * argv[]) {
    std::string igsFile = "data/intersection/case1-1.igs";
    try {
        NurbsFace face = nurbs_convert(igsFile);
        std::cout << "UBegin: " << face.get_ubegin() << std::endl;
        std::cout << "UEnd: " << face.get_uend() << std::endl;
        std::cout << "VBegin: " << face.get_vbegin() << std::endl;
        std::cout << "VEnd: " << face.get_vend() << std::endl;
        FILE * fid = fopen("../example/data/inter_case4-1_data.cu","w");
        if(fid == NULL)
        {
            printf("File open failed！\n");
            return 1;
        }
        fprintf(fid, "#include "nurbs_data.hpp"\n");
        fprintf(fid,"static int ivec[7]={");
        for (int i = 0; i < 7; i++) {
            if(i > 0){
                fprintf(fid,", ");
            }
            fprintf(fid,"%d", face.ivec[i]);
        }
        fprintf(fid,"};\n");
        fprintf(fid,"static int ndata = %d;\n", face.ndata);
        fprintf(fid,"static float data[%d]={", face.ndata);
        for (int i = 0; i < len; i++) {
            fprintf(fid, "%lf,", face.data[i]);
        }
        fprintf(fid,"};\n");
        fprintf(fid,"struct NurbsFace inter_case4_obj2 = {\n    .ivec = ivec,\n    .ndata = ndata,\n    .data = data\n};");

        fclose(fid);

        delete[] face.data;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}