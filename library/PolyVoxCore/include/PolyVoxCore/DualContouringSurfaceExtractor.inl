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
#include "SimpleVolume.h"

#include "PolyVoxCore/SurfaceMesh.h" 
#include "PolyVoxCore/VertexTypes.h" 

#include "Impl/QEF.h"

namespace PolyVox
{
	namespace
	{
		struct CellData
		{
			uint32_t vertexIndex;
			bool x;
			bool y;
			bool z;
			bool hasVertex;
		};
		
		template<typename VoxelType>
		struct EdgeData
		{
			VoxelType intersectionDistance;
			VoxelType edgeLength;
			bool intersects;
			Vector3DFloat normal;
		};
		
		template<typename SamplerType>
		Vector3DFloat computeCentralDifferenceGradient(const SamplerType& volIter)
		{
			//FIXME - Should actually use DensityType here, both in principle and because the maths may be
			//faster (and to reduce casts). So it would be good to add a way to get DensityType from a voxel.
			//But watch out for when the DensityType is unsigned and the difference could be negative.
			float voxel1nx = static_cast<float>(volIter.peekVoxel1nx0py0pz());
			float voxel1px = static_cast<float>(volIter.peekVoxel1px0py0pz());

			float voxel1ny = static_cast<float>(volIter.peekVoxel0px1ny0pz());
			float voxel1py = static_cast<float>(volIter.peekVoxel0px1py0pz());

			float voxel1nz = static_cast<float>(volIter.peekVoxel0px0py1nz());
			float voxel1pz = static_cast<float>(volIter.peekVoxel0px0py1pz());

			return Vector3DFloat
			(
				voxel1nx - voxel1px,
				voxel1ny - voxel1py,
				voxel1nz - voxel1pz
			);
		}
		
		template<typename VoxelType, typename ThresholdType>
		EdgeData<VoxelType> calculateEdge(VoxelType vA, VoxelType vB, Vector3DFloat gA, Vector3DFloat gB, ThresholdType threshold)
		{
			EdgeData<VoxelType> edge;
			edge.intersectionDistance = vA - threshold;
			edge.edgeLength = vA - vB;
			
			auto fraction = static_cast<float>(edge.intersectionDistance) / static_cast<float>(edge.edgeLength);
			if(fraction > 1.0 || fraction < 0.0)
			{
				edge.intersects = false;
				return edge;
			}
			else
			{
				edge.intersects = true;
			}
			
			edge.normal = ((gA * static_cast<float>(edge.edgeLength - edge.intersectionDistance) +  gB * static_cast<float>(edge.intersectionDistance)) / static_cast<float>(2 * edge.edgeLength));
			if(edge.normal.lengthSquared() > 0.000001f) 
			{
				edge.normal.normalise();
			}
			
			return edge;
		}
		
		template<typename VoxelType>
		PositionMaterialNormal computeVertex(std::array<EdgeData<VoxelType>, 12> edges)
		{
			Vector3DFloat massPoint{0,0,0}; //The average of the intersection vertices
			
			std::array<Vector3DFloat, 12> vertices;
			
			vertices[0] = {static_cast<float>(edges[0].intersectionDistance) / static_cast<float>(edges[0].edgeLength), 0, 0};
			vertices[1] = {0, static_cast<float>(edges[1].intersectionDistance) / static_cast<float>(edges[1].edgeLength), 0};
			vertices[2] = {0, 0, static_cast<float>(edges[2].intersectionDistance) / static_cast<float>(edges[2].edgeLength)};
			vertices[3] = {1, static_cast<float>(edges[3].intersectionDistance) / static_cast<float>(edges[3].edgeLength), 0};
			vertices[4] = {1, 0, static_cast<float>(edges[4].intersectionDistance) / static_cast<float>(edges[4].edgeLength)};
			vertices[5] = {0, 1, static_cast<float>(edges[5].intersectionDistance) / static_cast<float>(edges[5].edgeLength)};
			vertices[6] = {static_cast<float>(edges[6].intersectionDistance) / static_cast<float>(edges[6].edgeLength), 1, 0};
			vertices[7] = {static_cast<float>(edges[7].intersectionDistance) / static_cast<float>(edges[7].edgeLength), 0, 1};
			vertices[8] = {0, static_cast<float>(edges[8].intersectionDistance) / static_cast<float>(edges[8].edgeLength), 1};
			vertices[9] = {1, 1, static_cast<float>(edges[9].intersectionDistance) / static_cast<float>(edges[9].edgeLength)};
			vertices[10] = {1, static_cast<float>(edges[10].intersectionDistance) / static_cast<float>(edges[10].edgeLength), 1};
			vertices[11] = {static_cast<float>(edges[11].intersectionDistance) / static_cast<float>(edges[11].edgeLength), 1, 1};
			
			int numIntersections = 0;
			for(int i = 0; i < 12; ++i)
			{
				if(!edges[i].intersects)
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
				if(!edges[i].intersects)
				{
					continue;
				}
				
				Vector3DFloat normal = edges[i].normal;
				matrix[rows][0] = normal.getX();
				matrix[rows][1] = normal.getY();
				matrix[rows][2] = normal.getZ();
				
				Vector3DFloat p = vertices[i] - massPoint;
				Vector3DFloat product = normal * p;
				
				vector[rows] = product.getX() + product.getY() + product.getZ();
				
				cellVertexNormal += normal;
				
				++rows;
			}
			
			auto qefResult = evaluateQEF(matrix, vector, rows);
			Vector3DFloat vertexPosition = qefResult + massPoint;
			
			if(cellVertexNormal.lengthSquared() > 0.000001f) 
			{
				cellVertexNormal.normalise();
			}
			
			return {vertexPosition, cellVertexNormal, 0};
		}
	}
	
	template<typename VolumeType>
	SurfaceMesh<PositionMaterialNormal> dualContouringSurfaceExtractor(VolumeType* volData, Region region)
	{
		Timer timer;
		
		Region cellRegion{Vector3DInt32{0,0,0}, region.getDimensionsInVoxels()};
		SimpleVolume<CellData> cells{cellRegion}; //Cells should be 1 bigger than 'region' in each dimension
		
		typename VolumeType::Sampler volSampler{volData};
		volSampler.setPosition(region.getLowerCorner());
		volSampler.setWrapMode(WrapMode::Border, -100.0); // -100.0 is well below the threshold
		
		float threshold = 0;
		
		SurfaceMesh<PositionMaterialNormal> mesh;
		
		std::array<EdgeData<typename VolumeType::VoxelType>, 12> edges; //Create this now but it will be overwritten for each cell
		
		for(int32_t cellX = 0; cellX < cellRegion.getDimensionsInVoxels().getX(); cellX++)
		{
			for(int32_t cellY = 0; cellY < cellRegion.getDimensionsInVoxels().getY(); cellY++)
			{
				for(int32_t cellZ = 0; cellZ < cellRegion.getDimensionsInVoxels().getZ(); cellZ++)
				{
					//For each cell, calculate the vertex position
					volSampler.setPosition(region.getLowerCorner() + Vector3DInt32{cellX, cellY, cellZ} - Vector3DInt32{1,1,1});
					
					typename VolumeType::Sampler gradientSampler{volData};
					
					typename VolumeType::VoxelType v000{volSampler.getVoxel()};
					typename VolumeType::VoxelType v100{volSampler.peekVoxel1px0py0pz()};
					typename VolumeType::VoxelType v010{volSampler.peekVoxel0px1py0pz()};
					typename VolumeType::VoxelType v110{volSampler.peekVoxel1px1py0pz()};
					typename VolumeType::VoxelType v001{volSampler.peekVoxel0px0py1pz()};
					typename VolumeType::VoxelType v101{volSampler.peekVoxel1px0py1pz()};
					typename VolumeType::VoxelType v011{volSampler.peekVoxel0px1py1pz()};
					typename VolumeType::VoxelType v111{volSampler.peekVoxel1px1py1pz()};
					
					gradientSampler.setPosition(volSampler.getPosition());
					Vector3DFloat g000 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{1,0,0});
					Vector3DFloat g100 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{0,1,0});
					Vector3DFloat g010 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{1,1,0});
					Vector3DFloat g110 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{0,0,1});
					Vector3DFloat g001 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{1,0,1});
					Vector3DFloat g101 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{0,1,1});
					Vector3DFloat g011 = computeCentralDifferenceGradient(gradientSampler);
					gradientSampler.setPosition(volSampler.getPosition() + Vector3DInt32{1,1,1});
					Vector3DFloat g111 = computeCentralDifferenceGradient(gradientSampler);
					
					edges[0] = calculateEdge(v000, v100, g000, g100, threshold);
					edges[1] = calculateEdge(v000, v010, g000, g010, threshold);
					edges[2] = calculateEdge(v000, v001, g000, g001, threshold);
					
					edges[3] = calculateEdge(v100, v110, g100, g110, threshold);
					edges[4] = calculateEdge(v100, v101, g100, g101, threshold);
					
					edges[5] = calculateEdge(v010, v011, g010, g011, threshold);
					edges[6] = calculateEdge(v010, v110, g010, g110, threshold);
					
					edges[7] = calculateEdge(v001, v101, g001, g101, threshold);
					edges[8] = calculateEdge(v001, v011, g001, g011, threshold);
					
					edges[9] = calculateEdge(v110, v111, g110, g111, threshold);
					
					edges[10] = calculateEdge(v101, v111, g101, g111, threshold);
					
					edges[11] = calculateEdge(v011, v111, g011, g111, threshold);
					
					CellData cell;
					cell.x = edges[0].intersects;
					cell.y = edges[1].intersects;
					cell.z = edges[2].intersects;
					
					if(edges[0].intersects || edges[1].intersects || edges[2].intersects || edges[3].intersects || edges[4].intersects || edges[5].intersects || edges[6].intersects || edges[7].intersects || edges[8].intersects || edges[9].intersects || edges[10].intersects || edges[11].intersects)
					{
						cell.hasVertex = true;
						
						auto vertex = computeVertex(edges);
						
						vertex.setPosition(vertex.getPosition() + Vector3DFloat{static_cast<float>(cellX), static_cast<float>(cellY), static_cast<float>(cellZ)});
						
						mesh.addVertex(vertex);
						
						cell.vertexIndex = mesh.getNoOfVertices() - 1; //The index is of the last-added vertex
					}
					else
					{
						cell.hasVertex = false;
					}
					cells.setVoxel(cellX, cellY, cellZ, cell);
				}
			}
		}
		
		SimpleVolume<CellData>::Sampler cellSampler{&cells};
		cellSampler.setPosition(0,0,0);
		cellSampler.setWrapMode(WrapModes::Validate);
		
		for(int32_t cellX = 1; cellX < cellRegion.getDimensionsInVoxels().getX(); cellX++)
		{
			for(int32_t cellY = 1; cellY < cellRegion.getDimensionsInVoxels().getY(); cellY++)
			{
				for(int32_t cellZ = 1; cellZ < cellRegion.getDimensionsInVoxels().getZ(); cellZ++)
				{
					cellSampler.setPosition(cellX, cellY, cellZ);
					
					auto cell = cellSampler.getVoxel();
					
					if(cell.x)
					{
						auto v0 = cell;
						auto v1 = cellSampler.peekVoxel0px1ny0pz();
						auto v2 = cellSampler.peekVoxel0px0py1nz();
						auto v3 = cellSampler.peekVoxel0px1ny1nz();
						mesh.addTriangle(v0.vertexIndex, v1.vertexIndex, v2.vertexIndex);
						mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
					}
					
					if(cell.y)
					{
						auto v0 = cell;
						auto v1 = cellSampler.peekVoxel1nx0py0pz();
						auto v2 = cellSampler.peekVoxel0px0py1nz();
						auto v3 = cellSampler.peekVoxel1nx0py1nz();
						mesh.addTriangle(v0.vertexIndex, v1.vertexIndex, v2.vertexIndex);
						mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
					}
					
					if(cell.z)
					{
						auto v0 = cell;
						auto v1 = cellSampler.peekVoxel1nx0py0pz();
						auto v2 = cellSampler.peekVoxel0px1ny0pz();
						auto v3 = cellSampler.peekVoxel1nx1ny0pz();
						mesh.addTriangle(v0.vertexIndex, v1.vertexIndex, v2.vertexIndex);
						mesh.addTriangle(v3.vertexIndex, v2.vertexIndex, v1.vertexIndex);
					}
				}
			}
		}
		
		logTrace() << "Dual contouring surface extraction took " << timer.elapsedTimeInMilliSeconds() << "ms (Region size = " << region.getWidthInVoxels() << "x" << region.getHeightInVoxels() << "x" << region.getDepthInVoxels() << ")";
		
		return mesh;
	}
}
