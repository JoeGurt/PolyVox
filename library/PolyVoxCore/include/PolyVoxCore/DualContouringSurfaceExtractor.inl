/*******************************************************************************
Copyright (c) 2014 David Williams and Matt Williams

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
    distribution.
*******************************************************************************/
#include "PolyVoxCore/SurfaceMesh.h" 
#include "PolyVoxCore/VertexTypes.h" 

#include "Impl/Timer.h"

#include "Impl/QEF.h"

#include <type_traits>

//BUG We will get duplucation of edges if the surface is along region boundaries

namespace PolyVox
{
	namespace
	{
		template<typename VoxelType>
		struct EdgeData
		{
			EdgeData() : intersects(false) {}
			Vector3DFloat normal;
			float fraction; ///<fraction (0.0-1.0) along the edge in the positive direction that the intersection happens
			bool intersects;
		};
		
		template<typename VoxelType>
		struct CellData
		{
			EdgeData<VoxelType> edges[3];
			uint32_t vertexIndex;
		};
		
		template<typename VoxelType, typename ThresholdType>
		EdgeData<VoxelType> calculateEdge(const VoxelType& vA, const VoxelType& vB, const Vector3DFloat& gA, const Vector3DFloat& gB, const ThresholdType& threshold)
		{
			EdgeData<VoxelType> edge;
			
			edge.fraction = static_cast<float>(vA - threshold) / static_cast<float>(vA - vB);
			
			if(std::min(vA,vB) <= threshold && std::max(vA,vB) > threshold)
			{
				edge.intersects = true;
			}
			else
			{
				edge.intersects = false;
				return edge;
			}
			
			edge.normal = (gA * edge.fraction +  gB * (1.0f-edge.fraction));
			if(edge.normal.lengthSquared() != 0.0f) 
			{
				edge.normal.normalise();
			}
			
			return edge;
		}
		
		template<typename VoxelType>
		PositionMaterialNormal computeVertex(EdgeData<VoxelType>* edges[12])
		{
			Vector3DFloat massPoint{0,0,0}; //The average of the intersection vertices
			
			Vector3DFloat vertices[12];
			
			vertices[0] = {edges[0]->fraction, 0, 0};
			vertices[1] = {0, edges[1]->fraction, 0};
			vertices[2] = {0, 0, edges[2]->fraction};
			vertices[3] = {1, edges[3]->fraction, 0};
			vertices[4] = {1, 0, edges[4]->fraction};
			vertices[5] = {0, 1, edges[5]->fraction};
			vertices[6] = {edges[6]->fraction, 1, 0};
			vertices[7] = {edges[7]->fraction, 0, 1};
			vertices[8] = {0, edges[8]->fraction, 1};
			vertices[9] = {1, 1, edges[9]->fraction};
			vertices[10] = {1, edges[10]->fraction, 1};
			vertices[11] = {edges[11]->fraction, 1, 1};
			
			int numIntersections = 0;
			for(int i = 0; i < 12; ++i)
			{
				if(!edges[i]->intersects)
				{
					continue;
				}
				
				++numIntersections;
				massPoint += vertices[i];
			}
			
			massPoint /= numIntersections; //Make the average
			
			Vector3DFloat cellVertexNormal{0,0,0};
			
			double matrix[12][3];
			double vector[12];
			int rows = 0;
			
			for(int i = 0; i < 12; ++i)
			{
				if(!edges[i]->intersects)
				{
					continue;
				}
				
				Vector3DFloat normal = edges[i]->normal;
				matrix[rows][0] = normal.getX();
				matrix[rows][1] = normal.getY();
				matrix[rows][2] = normal.getZ();
				
				const Vector3DFloat p = vertices[i] - massPoint;
				const Vector3DFloat product = normal * p;
				
				vector[rows] = product.getX() + product.getY() + product.getZ();
				
				cellVertexNormal += normal;
				
				++rows;
			}
			
			const auto& vertexPosition = evaluateQEF(matrix, vector, rows) + massPoint;
			
			if(cellVertexNormal.lengthSquared() != 0.0f)
			{
				cellVertexNormal.normalise();
			}
			
			return {vertexPosition, cellVertexNormal, 0};
		}
		
		uint32_t convert(uint32_t x, uint32_t y, uint32_t z, uint32_t X, uint32_t Y)
		{
			return z*Y*X+y*X+x;
		}
	}
	
	template<typename VolumeType>
	SurfaceMesh<PositionMaterialNormal> dualContouringSurfaceExtractor(VolumeType* volData, Region region)
	{
		static_assert(std::is_signed<typename VolumeType::VoxelType>::value, "Voxel type must be signed");
		
		const float threshold = 0.0f;
		
		//Timer timer;
		Timer totalTimer;
		
		const auto regionXDimension = region.getDimensionsInVoxels().getX();
		const auto regionYDimension = region.getDimensionsInVoxels().getY();
		const auto regionZDimension = region.getDimensionsInVoxels().getZ();
		
		const auto gradientRegionXDimension = regionXDimension+2;
		const auto gradientRegionYDimension = regionYDimension+2;
		const auto gradientRegionZDimension = regionZDimension+2;
		
		std::vector<std::pair<const typename VolumeType::VoxelType, const Vector3DFloat>> gradients;
		gradients.reserve(gradientRegionXDimension * gradientRegionYDimension * gradientRegionZDimension);
		
		typename VolumeType::Sampler volSampler{volData};
		volSampler.setPosition(region.getLowerCorner() - Vector3DInt32{1,1,1});
		volSampler.setWrapMode(WrapMode::Border, -100.0); // -100.0 is well below the threshold
		
		const auto lowerCornerX = region.getLowerCorner().getX();
		const auto lowerCornerY = region.getLowerCorner().getZ();
		const auto lowerCornerZ = region.getLowerCorner().getX();
		
		//logTrace() << "Setup took " << timer.elapsedTimeInMilliSeconds();
		//timer.start();
		
		for(int32_t z = 0; z < gradientRegionZDimension; z++)
		{
			volSampler.setPosition(lowerCornerX-1, lowerCornerY-1, lowerCornerZ+z-1); //Reset x and y and increment z
			for(int32_t y = 0; y < gradientRegionYDimension; y++)
			{
				volSampler.setPosition(lowerCornerX-1, lowerCornerY+y-1, lowerCornerZ+z-1); //Reset x and increment y (z remains the same)
				for(int32_t x = 0; x < gradientRegionXDimension; x++)
				{
					volSampler.movePositiveX(); //Increment x
					
					const auto& voxel = volSampler.getVoxel();
					const auto& voxel1px = volSampler.peekVoxel1px0py0pz();
					const auto& voxel1py = volSampler.peekVoxel0px1py0pz();
					const auto& voxel1pz = volSampler.peekVoxel0px0py1pz();
					
					const auto& voxel1nx = volSampler.peekVoxel1nx0py0pz();
					const auto& voxel1ny = volSampler.peekVoxel0px1ny0pz();
					const auto& voxel1nz = volSampler.peekVoxel0px0py1nz();
					
					gradients.emplace_back(voxel, Vector3DFloat(voxel1nx - voxel1px, voxel1ny - voxel1py, voxel1nz - voxel1pz));
				}
			}
		}
		
		//logTrace() << "Gradients took " << timer.elapsedTimeInMilliSeconds();
		//timer.start();
		
		const auto cellRegionXDimension = regionXDimension+2;
		const auto cellRegionYDimension = regionYDimension+2;
		const auto cellRegionZDimension = regionZDimension+2;
		
		std::vector<CellData<typename VolumeType::VoxelType>> cells;
		cells.reserve(cellRegionXDimension * cellRegionYDimension * cellRegionZDimension);
		
		for(int32_t cellZ = 0; cellZ < cellRegionZDimension; cellZ++)
		{
			for(int32_t cellY = 0; cellY < cellRegionYDimension; cellY++)
			{
				for(int32_t cellX = 0; cellX < cellRegionXDimension; cellX++)
				{
					//For each cell, calculate the edge intersection points and normals
					const auto& g000 = gradients[convert(cellX, cellY, cellZ, cellRegionXDimension, cellRegionYDimension)];
					
					//For the last columns/rows, only calculate the interior edge
					if(cellX < cellRegionXDimension-1 && cellY < cellRegionYDimension-1 && cellZ < cellRegionZDimension-1) //This is the main bulk
					{
						const auto& g100 = gradients[convert(cellX+1, cellY, cellZ, cellRegionXDimension, cellRegionYDimension)];
						const auto& g010 = gradients[convert(cellX, cellY+1, cellZ, cellRegionXDimension, cellRegionYDimension)];
						const auto& g001 = gradients[convert(cellX, cellY, cellZ+1, cellRegionXDimension, cellRegionYDimension)];
						cells.push_back({calculateEdge(g000.first, g100.first, g000.second, g100.second, threshold), calculateEdge(g000.first, g010.first, g000.second, g010.second, threshold), calculateEdge(g000.first, g001.first, g000.second, g001.second, threshold)});
					}
					else if(cellX == cellRegionXDimension-1 || cellY == cellRegionYDimension-1 || cellZ == cellRegionZDimension-1) //This is the three far edges and the far corner
					{
						cells.push_back({}); //Default and empty
					}
					else if(cellX == cellRegionXDimension-1) //Far x side
					{
						const auto& g100 = gradients[convert(cellX+1, cellY, cellZ, cellRegionXDimension, cellRegionYDimension)];
						cells.push_back({calculateEdge(g000.first, g100.first, g000.second, g100.second, threshold), EdgeData<typename VolumeType::VoxelType>(), EdgeData<typename VolumeType::VoxelType>()});
					}
					else if(cellY == cellRegionYDimension-1) //Far y side
					{
						const auto& g010 = gradients[convert(cellX+1, cellY, cellZ, cellRegionXDimension, cellRegionYDimension)];
						cells.push_back({EdgeData<typename VolumeType::VoxelType>(), calculateEdge(g000.first, g010.first, g000.second, g010.second, threshold), EdgeData<typename VolumeType::VoxelType>()});
					}
					else if(cellZ == cellRegionZDimension-1) //Far z side
					{
						const auto& g001 = gradients[convert(cellX+1, cellY, cellZ, cellRegionXDimension, cellRegionYDimension)];
						cells.push_back({EdgeData<typename VolumeType::VoxelType>(), EdgeData<typename VolumeType::VoxelType>(), calculateEdge(g000.first, g001.first, g000.second, g001.second, threshold)});
					}
				}
			}
		}
		
		//logTrace() << "Edges took " << timer.elapsedTimeInMilliSeconds();
		//timer.start();
		
		EdgeData<typename VolumeType::VoxelType>* edges[12]; //Create this now but it will be overwritten for each cell
		
		SurfaceMesh<PositionMaterialNormal> mesh;
		
		for(int32_t cellZ = 0; cellZ < cellRegionZDimension; cellZ++)
		{
			for(int32_t cellY = 0; cellY < cellRegionYDimension; cellY++)
			{
				for(int32_t cellX = 0; cellX < cellRegionXDimension; cellX++)
				{
					if(cellZ >= 1 && cellY >= 1 && cellX >= 1)
					{
						//After the first rows and columns are done, start calculating vertex positions
						const int32_t cellXVertex = cellX-1;
						const int32_t cellYVertex = cellY-1;
						const int32_t cellZVertex = cellZ-1;
						
						auto& cell = cells[convert(cellXVertex, cellYVertex, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
						
						edges[0] = &cell.edges[0];
						edges[1] = &cell.edges[1];
						edges[2] = &cell.edges[2];
						
						edges[3] = &cells[convert(cellXVertex+1, cellYVertex, cellZVertex, cellRegionXDimension, cellRegionYDimension)].edges[1];
						edges[4] = &cells[convert(cellXVertex+1, cellYVertex, cellZVertex, cellRegionXDimension, cellRegionYDimension)].edges[2];
						
						edges[5] = &cells[convert(cellXVertex, cellYVertex+1, cellZVertex, cellRegionXDimension, cellRegionYDimension)].edges[2];
						edges[6] = &cells[convert(cellXVertex, cellYVertex+1, cellZVertex, cellRegionXDimension, cellRegionYDimension)].edges[0];
						
						edges[7] = &cells[convert(cellXVertex, cellYVertex, cellZVertex+1, cellRegionXDimension, cellRegionYDimension)].edges[0];
						edges[8] = &cells[convert(cellXVertex, cellYVertex, cellZVertex+1, cellRegionXDimension, cellRegionYDimension)].edges[1];
						
						edges[9] = &cells[convert(cellXVertex+1, cellYVertex+1, cellZVertex, cellRegionXDimension, cellRegionYDimension)].edges[2];
						
						edges[10] = &cells[convert(cellXVertex+1, cellYVertex, cellZVertex+1, cellRegionXDimension, cellRegionYDimension)].edges[1];
						
						edges[11] = &cells[convert(cellXVertex, cellYVertex+1, cellZVertex+1, cellRegionXDimension, cellRegionYDimension)].edges[0];
						
						if(edges[0]->intersects || edges[1]->intersects || edges[2]->intersects || edges[3]->intersects || edges[4]->intersects || edges[5]->intersects || edges[6]->intersects || edges[7]->intersects || edges[8]->intersects || edges[9]->intersects || edges[10]->intersects || edges[11]->intersects) //'if' Maybe not needed?
						{
							auto vertex = computeVertex(edges);
							
							vertex.setPosition({vertex.getPosition().getX()+cellXVertex, vertex.getPosition().getY()+cellYVertex, vertex.getPosition().getZ()+cellZVertex});
							
							cell.vertexIndex = mesh.addVertex(vertex);
							
							if(cellZVertex >= 1 && cellYVertex >= 1 && cellXVertex >= 1)
							{
								//Once the second rows and colums are done, start connecting up edges
								if(cell.edges[0].intersects)
								{
									const auto& v1 = cells[convert(cellXVertex, cellYVertex-1, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
									const auto& v2 = cells[convert(cellXVertex, cellYVertex, cellZVertex-1, cellRegionXDimension, cellRegionYDimension)];
									const auto& v3 = cells[convert(cellXVertex, cellYVertex-1, cellZVertex-1, cellRegionXDimension, cellRegionYDimension)];
									mesh.addTriangle(cell.vertexIndex, v1.vertexIndex, v2.vertexIndex);
									mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
								}
								
								if(cell.edges[1].intersects)
								{
									const auto& v1 = cells[convert(cellXVertex-1, cellYVertex, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
									const auto& v2 = cells[convert(cellXVertex, cellYVertex, cellZVertex-1, cellRegionXDimension, cellRegionYDimension)];
									const auto& v3 = cells[convert(cellXVertex-1, cellYVertex, cellZVertex-1, cellRegionXDimension, cellRegionYDimension)];
									mesh.addTriangle(cell.vertexIndex, v1.vertexIndex, v2.vertexIndex);
									mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
								}
								
								if(cell.edges[2].intersects)
								{
									const auto& v1 = cells[convert(cellXVertex-1, cellYVertex, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
									const auto& v2 = cells[convert(cellXVertex, cellYVertex-1, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
									const auto& v3 = cells[convert(cellXVertex-1, cellYVertex-1, cellZVertex, cellRegionXDimension, cellRegionYDimension)];
									mesh.addTriangle(cell.vertexIndex, v1.vertexIndex, v2.vertexIndex);
									mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
								}
							}
						}
					}
				}
			}
		}
		
		//logTrace() << "Vertices and quads took " << timer.elapsedTimeInMilliSeconds();
		//timer.start();
		
		logTrace() << "Dual contouring surface extraction took " << totalTimer.elapsedTimeInMilliSeconds() << "ms (Region size = " << region.getWidthInVoxels() << "x" << region.getHeightInVoxels() << "x" << region.getDepthInVoxels() << ")";
		
		logTrace() << mesh.getNoOfVertices();
		
		return mesh;
	}
}
