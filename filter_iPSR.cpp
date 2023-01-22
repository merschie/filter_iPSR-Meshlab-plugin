/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History
$Log: samplefilter.cpp,v $
Revision 1.3  2006/11/29 00:59:20  cignoni
Cleaned plugins interface; changed useless help class into a plain string

Revision 1.2  2006/11/27 06:57:21  cignoni
Wrong way of using the __DATE__ preprocessor symbol

Revision 1.1  2006/09/25 09:24:39  e_cerisoli
add samplefilter

****************************************************************************/

// #include <math.h>
// #include <stdlib.h>
// #include <time.h>

// #include <vcg/complex/algorithms/create/platonic.h>

#include "filter_iPSR.h"
// #include "src/Geometry.h"
// #include "src/PoissonParam.h"



#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include "src/kdtree.h"
#include "src/utility.h"
#include "src/PoissonRecon.h"
#include "src/PointStreamData.h"

using namespace std;
using namespace vcg;

//int Execute2(PoissonParam &Par, vector<Point3D<float> > Pts, vector<Point3D<float> > Nor, 	CoredVectorMeshData &mesh, Point3D<float> &newCenter, float &newScale, vcg::CallBackPos *cb );

// Constructor usually performs only two simple tasks of filling the two lists
//  - typeList: with all the possible id of the filtering actions
//  - actionList with the corresponding actions. If you want to add icons to your filtering actions you can do here by construction the QActions accordingly

PoissonPlugin::PoissonPlugin()
{
	typeList = {FP_POISSON_RECON};

	foreach(ActionIDType tt , types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QString PoissonPlugin::pluginName() const
{
	return "FilteriPSR";
}

// ST() must return the very short string describing each filtering action
// (this string is used also to define the menu entry)
 QString PoissonPlugin::filterName(ActionIDType filterId) const
{
  switch(filterId) {
  case FP_POISSON_RECON :  return QString("Iterative Poisson Surface Reconstruction");
        default : assert(0);
    }
    return QString("Error: Unknown Filter");
}

// Info() must return the longer string describing each filtering action
// (this string is used in the About plugin dialog)
 QString PoissonPlugin::filterInfo(ActionIDType filterId) const
{
  switch(filterId) {
        case FP_POISSON_RECON :  return QString("Iterative Poisson Surface Reconstruction (iPSR) for Unoriented Points ");
        default : assert(0);
    }
    return QString("Error: Unknown Filter");
}

// This function define the needed parameters for each filter. Return true if the filter has some parameters
// it is called every time, so you can set the default value of parameters according to the mesh
// For each parmeter you need to define,
// - the name of the parameter,
// - the string shown in the dialog
// - the default value
// - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
RichParameterList PoissonPlugin::initParameterList(const QAction* action, const MeshModel &)
{
     RichParameterList parlst;
     switch(ID(action))	 {
        case FP_POISSON_RECON :
          //parlst.addParam(new RichBool ("RecomputeNormals",
            //								false,
            //								"Recompute normals",
            //								"Do not use the current normals, but recompute them from scratch considering the vertices of the mesh as an unstructured point cloud.");
          //parlst.addParam(new RichBool ("UseConfidence",
            //								true,
            //								"Use Quality",
            //								"Use the per vertex quality as a confidence value\n");
            parlst.addParam(RichInt ("iters",
                                            30,
                                            "maximum number of iterations",
                                            "The maximum number of iterations. The default value of this parameter is 30.\n"));
            parlst.addParam(RichInt ("pointWeight",
                                            10,
                                            "interpolation weight",
                                            "The pointWeight parameter of screened Poisson surface reconstruction. The default value for this parameter is 10.\n"));
            parlst.addParam(RichInt ("depth",
                                            10,
                                            "Samples per Node",
                                            "The depth parameter of screened Poisson surface reconstruction. It is the maximum possible depth of the octree. The default value of this parameter is 10.\n"));
            parlst.addParam(RichFloat ("neighbors",
                                             10,
                                             "number of neighbors",
                                             "The number of the closest sample points to search from every reconstructed triangle face. The suggested range is between 10 and 20. The default value of this parameter is 10.\n"));

            break;
   default: break; // do not add any parameter for the other filters
  }
     return parlst;
}

// The Real Core Function doing the actual mesh processing.
// Move Vertex of a random quantity
std::map<std::string, QVariant> PoissonPlugin::applyFilter(
		const QAction * action,
		const RichParameterList & par,
		MeshDocument &md,
		unsigned int&,
		vcg::CallBackPos *cb)
{
	typedef double REAL;
	const unsigned int DIM = 3U;
    vector<pair<Point<double, 3>, NormaliSPR<double, 3>>> points_normals;
	//cancel read in, use mesh from meshlab
	//ply_reader<REAL, DIM>(input_name, points_normals);
	
	//parse the parameters
	int iters =par.getInt("iters");
	int pointweight =par.getInt("pointWeight");
	int depth =par.getInt("depth");
	int k_neighbors =par.getInt("neighbors");


	string command = "PoissonRecon --in i.ply --out o.ply --bType 2 --depth " + to_string(depth) + " --pointWeight " + to_string(pointweight);
	log("command: %s",command.c_str());
	vector<string> cmd = split(command);
	vector<char *> argv_str(cmd.size());
	for (size_t i = 0; i < cmd.size(); ++i)
		argv_str[i] = &cmd[i][0];

	XForm<REAL, DIM + 1> iXForm;
	vector<double> weight_samples;
	// sample points by the octree
	points_normals = sample_points<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, iXForm, &weight_samples);

	// initialize normals randomly
	printf("random initialization...\n");
	NormaliSPR<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
	srand(0);
	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		do
		{
			points_normals[i].second = Point<REAL, DIM>(rand() % 1001 - 500.0, rand() % 1001 - 500.0, rand() % 1001 - 500.0);
		} while (points_normals[i].second == zero_normal);
		normalize<REAL, DIM>(points_normals[i].second);
	}

	// construct the Kd-Tree
	kdt::KDTree<kdt::KDTreePoint> tree;
	{
		vector<kdt::KDTreePoint> vertices;
		vertices.reserve(points_normals.size());
		for (size_t i = 0; i < points_normals.size(); ++i)
		{
			array<double, 3> p{points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]};
			vertices.push_back(kdt::KDTreePoint(p));
		}
		tree.build(vertices);
	}

	pair<vector<Point<REAL, DIM>>, vector<vector<int>>> mesh;

	// iterations
	int epoch = 0;
	while (epoch < iters)
	{
		++epoch;
		printf("Iter: %d\n", epoch);

		vector<Point<REAL, DIM>>().swap(mesh.first);
		vector<vector<int>>().swap(mesh.second);

		// Poisson reconstruction
		mesh = poisson_reconstruction<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, &weight_samples);

		vector<vector<int>> nearestSamples(mesh.second.size());
		vector<Point<REAL, DIM>> normals(mesh.second.size());

		// compute face normals and map them to sample points
#pragma omp parallel for
		for (int i = 0; i < (int)nearestSamples.size(); i++)
		{
			if (mesh.second[i].size() == 3)
			{
				Point<REAL, DIM> c = mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]];
				c /= 3;
				array<REAL, DIM> a{c[0], c[1], c[2]};
				nearestSamples[i] = tree.knnSearch(kdt::KDTreePoint(a), k_neighbors);
				normals[i] = Point<REAL, DIM>::CrossProduct(mesh.first[mesh.second[i][1]] - mesh.first[mesh.second[i][0]], mesh.first[mesh.second[i][2]] - mesh.first[mesh.second[i][0]]);
			}
		}

		// update sample point normals
		vector<NormaliSPR<REAL, DIM>> projective_normals(points_normals.size(), zero_normal);
		for (size_t i = 0; i < nearestSamples.size(); i++)
		{
			for (size_t j = 0; j < nearestSamples[i].size(); ++j)
			{
				projective_normals[nearestSamples[i][j]].normal[0] += normals[i][0];
				projective_normals[nearestSamples[i][j]].normal[1] += normals[i][1];
				projective_normals[nearestSamples[i][j]].normal[2] += normals[i][2];
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (int)projective_normals.size(); ++i)
			normalize<REAL, DIM>(projective_normals[i]);

		// compute the average normal variation of the top 1/1000 points
		size_t heap_size = static_cast<size_t>(ceil(points_normals.size() / 1000.0));
		priority_queue<double, vector<double>, greater<double>> min_heap;
		for (size_t i = 0; i < points_normals.size(); ++i)
		{
			if (!(projective_normals[i] == zero_normal))
			{
				double diff = Point<REAL, DIM>::SquareNorm((projective_normals[i] - points_normals[i].second).normal);
				if (min_heap.size() < heap_size)
					min_heap.push(diff);
				else if (diff > min_heap.top())
				{
					min_heap.pop();
					min_heap.push(diff);
				}

				points_normals[i].second = projective_normals[i];
			}
		}

		heap_size = min_heap.size();
		double ave_max_diff = 0;
		while (!min_heap.empty())
		{
			ave_max_diff += sqrt(min_heap.top());
			min_heap.pop();
		}
		ave_max_diff /= heap_size;
		printf("normals variation %f\n", ave_max_diff);
		if (ave_max_diff < 0.175)
			break;
	}

	mesh = poisson_reconstruction<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, &weight_samples);

	//output_ply(output_name, mesh, iXForm);



	
	// output_sample_points_and_normals<REAL, DIM>("points_normals_samples.ply", points_normals, iXForm);
	// output_all_points_and_normals<REAL, DIM>("points_normals_all.ply", input_name, points_normals, tree, iXForm);





	// if (ID(action) == FP_POISSON_RECON) {
	// 	MeshModel &m=*md.mm();
	// 	MeshModel &pm =*md.addNewMesh("","Poisson mesh");
	// 	vector<Point3D<float> > Pts(m.cm.vn);
	// 	vector<Point3D<float> > Nor(m.cm.vn);
	// 	CoredVectorMeshData mesh;
	
	// 	if (m.hasDataMask(MeshModel::MM_WEDGTEXCOORD)){
	// 		m.clearDataMask(MeshModel::MM_WEDGTEXCOORD);
	// 	}
	// 	if (m.hasDataMask(MeshModel::MM_VERTTEXCOORD)){
	// 		m.clearDataMask(MeshModel::MM_VERTTEXCOORD);
	// 	}
	
	// 	//Useless control on the normals. It can just avoid crashes derived from an importer setting up to [0.0f,0.0f,0.0f] the normal vectors of a mesh without per-vertex normal attribute.
	// 	int zeronrm = 0;
	// 	for(CMeshO::VertexIterator vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi)
	// 	{
	// 		if(!(*vi).IsD())
	// 		{
	// 			if ((*vi).N() == Point3m(0.0f,0.0f,0.0f))
	// 				++zeronrm;
	// 		}
	// 	}
	
	// 	if (zeronrm == m.cm.vn)
	// 	{
	// 		log(GLLogStream::SYSTEM,"All the normal vectors are set to [0.0,0.0,0.0]. Poisson reconstruction filter requires a set of valid per-vertex normal. Filter will be aborted.");
	// 		throw MLException("All the normal vectors are set to [0.0,0.0,0.0]. Poisson reconstruction filter requires a set of valid per-vertex normal. Filter will be aborted.");
	// 	}
	
	// 	int cnt=0;
	// 	for(CMeshO::VertexIterator vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi)
	// 	if(!(*vi).IsD()){
	// 			(*vi).N().Normalize();
	// 			for(int ii=0;ii<3;++ii){
	// 					Pts[cnt].coords[ii]=(*vi).P()[ii];
	// 					Nor[cnt].coords[ii]=(*vi).N()[ii];
	// 			}
	// 			++cnt;
	// 		}
	// 	assert(cnt==m.cm.vn);
	// 	// Log function dump textual info in the lower part of the MeshLab screen.
	// 	PoissonParam pp;
	// 	pp.Depth=par.getInt("OctDepth");
	// 	pp.SamplesPerNode = par.getFloat("SamplesPerNode");
	// 	pp.SolverDivide=par.getInt("SolverDivide");
	// 	pp.Offset = par.getFloat("Offset");
	// 	Point3D<float> center;
	// 	float scale;
	
	// 	Execute2(pp, Pts, Nor, mesh,center,scale,cb);
	// 	mesh.resetIterator();
	// 	int vm = mesh.outOfCorePointCount()+mesh.inCorePoints.size();
	// 	int fm = mesh.triangleCount();
	
	// 	log("Successfully created a mesh of %i vert and %i faces",vm,fm);
	
	// 	//m.cm.Clear();
	
	// 	tri::Allocator<CMeshO>::AddVertices(pm.cm,vm);
	// 	tri::Allocator<CMeshO>::AddFaces(pm.cm,fm);
	
	//   Point3D<float> p;
	// 	int i;
	// 	for (i=0; i < int(mesh.inCorePoints.size()); i++){
	// 		p=mesh.inCorePoints[i];
	// 		pm.cm.vert[i].P()[0] = p.coords[0]*scale+center.coords[0];
	// 		pm.cm.vert[i].P()[1] = p.coords[1]*scale+center.coords[1];
	// 		pm.cm.vert[i].P()[2] = p.coords[2]*scale+center.coords[2];
	// 		}
	// 	for (int ii=0; ii < mesh.outOfCorePointCount(); ii++){
	// 		mesh.nextOutOfCorePoint(p);
	// 		pm.cm.vert[ii+i].P()[0] = p.coords[0]*scale+center.coords[0];
	// 		pm.cm.vert[ii+i].P()[1] = p.coords[1]*scale+center.coords[1];
	// 		pm.cm.vert[ii+i].P()[2] = p.coords[2]*scale+center.coords[2];
	// 	}
	
	// TriangleIndex tIndex;
	// int inCoreFlag;
	// int nr_faces=mesh.triangleCount();
	// for (i=0; i < nr_faces; i++){
	// 		//
	// 		// create and fill a struct that the ply code can handle
	// 		//
	// 		mesh.nextTriangle(tIndex,inCoreFlag);
	// 		if(!(inCoreFlag & CoredMeshData::IN_CORE_FLAG[0])){tIndex.idx[0]+=int(mesh.inCorePoints.size());}
	// 		if(!(inCoreFlag & CoredMeshData::IN_CORE_FLAG[1])){tIndex.idx[1]+=int(mesh.inCorePoints.size());}
	// 		if(!(inCoreFlag & CoredMeshData::IN_CORE_FLAG[2])){tIndex.idx[2]+=int(mesh.inCorePoints.size());}
	// 		for(int j=0; j < 3; j++)
	// 		{
	// 			pm.cm.face[i].V(j) = &pm.cm.vert[tIndex.idx[j]];
	// 		}
	// 		//ply_put_element(ply, (void *) &ply_face);
	// 		//delete[] ply_face.vertices;
	// 	}  // for, write faces
	
	
	// //	for(int i=0;i<mesh.inCorePoints.size();++i){
	// //		mesh.triangles[i].idx[0]+=mesh.inCorePoints.size();
	// //		mesh.triangles[i].idx[1]+=mesh.inCorePoints.size();
	// //		mesh.triangles[i].idx[2]+=mesh.inCorePoints.size();
	// //		}
	// //	Build(m.cm,mesh.inCorePoints,mesh.triangles);
	// 	log("Successfully created a mesh of %i faces",pm.cm.vn);
	
	// 	pm.updateBoxAndNormals();
	// }
	// else {
	// 	wrongActionCalled(action);
	// }
	return std::map<std::string, QVariant>();
}
 PoissonPlugin::FilterClass PoissonPlugin::getClass(const QAction *action) const
{
  switch(ID(action))
  {
    case FP_POISSON_RECON :
            return FilterClass (FilterPlugin::PointSet + FilterPlugin::Remeshing) ;
    default: assert(0);
  }
  return FilterClass(0);
}


MESHLAB_PLUGIN_NAME_EXPORTER(PoissonPlugin)
