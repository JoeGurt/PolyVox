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

//BUG We will get duplucation of edges if the surface is along region boundaries

namespace PolyVox
{
	namespace
	{
		template<typename VoxelType>
		struct EdgeData
		{
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
			
			if(edge.fraction > 1.0 || edge.fraction < 0.0)
			{
				edge.intersects = false;
				return edge;
			}
			else
			{
				edge.intersects = true;
			}
			
			edge.normal = (gA * edge.fraction +  gB * (1.0f-edge.fraction));
			if(edge.normal.lengthSquared() > 0.000001f) 
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
				
				Vector3DFloat p = vertices[i] - massPoint;
				const Vector3DFloat product = normal * p;
				
				vector[rows] = product.getX() + product.getY() + product.getZ();
				
				cellVertexNormal += normal;
				
				++rows;
			}
			
			const auto vertexPosition = evaluateQEF(matrix, vector, rows) + massPoint;
			
			if(cellVertexNormal.lengthSquared() > 0.000001f) 
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
		Timer timer;
		
		const auto regionXDimension = region.getDimensionsInVoxels().getX();
		const auto regionYDimension = region.getDimensionsInVoxels().getY();
		const auto regionZDimension = region.getDimensionsInVoxels().getZ();
		
		const auto cellRegionXDimension = regionXDimension+1;
		const auto cellRegionYDimension = regionYDimension+1;
		const auto cellRegionZDimension = regionZDimension+1;
		
		const auto gradientRegionXDimension = regionXDimension+2;
		const auto gradientRegionYDimension = regionYDimension+2;
		const auto gradientRegionZDimension = regionZDimension+2;
		
		std::vector<std::pair<typename VolumeType::VoxelType*, Vector3DFloat>> gradients;
		gradients.reserve(gradientRegionXDimension * gradientRegionYDimension * gradientRegionZDimension);
		
		std::vector<CellData<typename VolumeType::VoxelType>> cells;
		cells.reserve(cellRegionXDimension * cellRegionYDimension * cellRegionZDimension);
		
		typename VolumeType::Sampler volSampler{volData};
		volSampler.setPosition(region.getLowerCorner());
		volSampler.setWrapMode(WrapMode::Border, -100.0); // -100.0 is well below the threshold
		
		const float threshold = 0.0f;
		
		SurfaceMesh<PositionMaterialNormal> mesh;
		
		EdgeData<typename VolumeType::VoxelType>* edges[12]; //Create this now but it will be overwritten for each cell
		
		const auto lowerCornerX = region.getLowerCorner().getX();
		const auto lowerCornerY = region.getLowerCorner().getZ();
		const auto lowerCornerZ = region.getLowerCorner().getX();
		
		for(int32_t cellZ = 0; cellZ < cellRegionZDimension; cellZ++)
		{
			for(int32_t cellY = 0; cellY < cellRegionYDimension; cellY++)
			{
				for(int32_t cellX = 0; cellX < cellRegionXDimension; cellX++)
				{
					//For each cell, calculate the edge intersection points and normals
					volSampler.setPosition(lowerCornerX+cellX-1, lowerCornerY+cellY-1, lowerCornerZ+cellZ-1);
					
					const auto& voxel = static_cast<float>(volSampler.getVoxel());
					const auto& voxel1px = static_cast<float>(volSampler.peekVoxel1px0py0pz());
					const auto& voxel1py = static_cast<float>(volSampler.peekVoxel0px1py0pz());
					const auto& voxel1pz = static_cast<float>(volSampler.peekVoxel0px0py1pz());
					
					const auto& voxel1nx = static_cast<float>(volSampler.peekVoxel1nx0py0pz());
					const auto& voxel1ny = static_cast<float>(volSampler.peekVoxel0px1ny0pz());
					const auto& voxel1nz = static_cast<float>(volSampler.peekVoxel0px0py1nz());
					const Vector3DFloat g000(voxel1nx - voxel1px, voxel1ny - voxel1py, voxel1nz - voxel1pz);
					volSampler.movePositiveX();
					const auto& voxel2px = static_cast<float>(volSampler.peekVoxel1px0py0pz());
					const auto& voxel1px1ny = static_cast<float>(volSampler.peekVoxel0px1ny0pz());
					const auto& voxel1px1py = static_cast<float>(volSampler.peekVoxel0px1py0pz());
					const auto& voxel1px1nz = static_cast<float>(volSampler.peekVoxel0px0py1nz());
					const auto& voxel1px1pz = static_cast<float>(volSampler.peekVoxel0px0py1pz());
					const Vector3DFloat g100(voxel - voxel2px, voxel1px1ny - voxel1px1py, voxel1px1nz - voxel1px1pz);
					volSampler.moveNegativeX();
					volSampler.movePositiveY();
					const auto& voxel1nx1py = static_cast<float>(volSampler.peekVoxel1nx0py0pz());
					const auto& voxel2py = static_cast<float>(volSampler.peekVoxel0px1py0pz());
					const auto& voxel1py1nz = static_cast<float>(volSampler.peekVoxel0px0py1nz());
					const auto& voxel1py1pz = static_cast<float>(volSampler.peekVoxel0px0py1pz());
					const Vector3DFloat g010(voxel1nx1py - voxel1px1py, voxel - voxel2py, voxel1py1nz - voxel1py1pz);
					volSampler.moveNegativeY();
					volSampler.movePositiveZ();
					const auto& voxel1nx1pz = static_cast<float>(volSampler.peekVoxel1nx0py0pz());
					const auto& voxel1ny1pz = static_cast<float>(volSampler.peekVoxel0px1ny0pz());
					const auto& voxel2pz = static_cast<float>(volSampler.peekVoxel0px0py1pz());
					const Vector3DFloat g001(voxel1nx1pz - voxel1px1pz, voxel1ny1pz - voxel1py1pz, voxel - voxel2pz);
					
					cells.push_back({calculateEdge(voxel, voxel1px, g000, g100, threshold), calculateEdge(voxel, voxel1py, g000, g010, threshold), calculateEdge(voxel, voxel1pz, g000, g001, threshold)});
					
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
		
		//logTrace() << "Dual contouring surface extraction took " << timer.elapsedTimeInMilliSeconds() << "ms (Region size = " << region.getWidthInVoxels() << "x" << region.getHeightInVoxels() << "x" << region.getDepthInVoxels() << ")";
		
		logTrace() << "Dual contouring surface extraction took " << timer.elapsedTimeInMilliSeconds() << "ms (Region size = " << region.getWidthInVoxels() << "x" << region.getHeightInVoxels() << "x" << region.getDepthInVoxels() << ")";
		
		std::cout << mesh.getNoOfVertices() << std::endl;
		
		return mesh;
	}
}
