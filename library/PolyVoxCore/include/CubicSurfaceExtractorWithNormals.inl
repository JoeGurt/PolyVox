/*******************************************************************************
Copyright (c) 2005-2009 David Williams

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

#include "Array.h"
#include "MaterialDensityPair.h"
#include "SurfaceMesh.h"
#include "PolyVoxImpl/MarchingCubesTables.h"
#include "VertexTypes.h"

namespace PolyVox
{
	template <typename VoxelType>
	CubicSurfaceExtractorWithNormals<VoxelType>::CubicSurfaceExtractorWithNormals(Volume<VoxelType>* volData, Region region, SurfaceMesh<PositionMaterialNormal>* result)
		:m_volData(volData)
		,m_sampVolume(volData)
		,m_regSizeInVoxels(region)
		,m_meshCurrent(result)
	{
		m_regSizeInCells = m_regSizeInVoxels;
		m_regSizeInCells.setUpperCorner(m_regSizeInCells.getUpperCorner() - Vector3DInt16(1,1,1));

		m_meshCurrent->clear();
	}

	template <typename VoxelType>
	void CubicSurfaceExtractorWithNormals<VoxelType>::execute()
	{		
		for(int16_t z = m_regSizeInVoxels.getLowerCorner().getZ(); z < m_regSizeInVoxels.getUpperCorner().getZ(); z++)
		{
			for(int16_t y = m_regSizeInVoxels.getLowerCorner().getY(); y < m_regSizeInVoxels.getUpperCorner().getY(); y++)
			{
				for(int16_t x = m_regSizeInVoxels.getLowerCorner().getX(); x < m_regSizeInVoxels.getUpperCorner().getX(); x++)
				{
					uint16_t regX = x - m_regSizeInVoxels.getLowerCorner().getX();
					uint16_t regY = y - m_regSizeInVoxels.getLowerCorner().getY();
					uint16_t regZ = z - m_regSizeInVoxels.getLowerCorner().getZ();

					int currentVoxel = m_volData->getVoxelAt(x,y,z).getDensity() >= VoxelType::getThreshold();

					int plusXVoxel = m_volData->getVoxelAt(x+1,y,z).getDensity()  >= VoxelType::getThreshold();
					if(currentVoxel > plusXVoxel)
					{
						int material = m_volData->getVoxelAt(x,y,z).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ - 0.5f), Vector3DFloat(1.0f, 0.0f, 0.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(1.0f, 0.0f, 0.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(1.0f, 0.0f, 0.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(1.0f, 0.0f, 0.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v2,v1);
						m_meshCurrent->addTriangleCubic(v1,v2,v3);
					}
					if(currentVoxel < plusXVoxel)
					{
						int material = m_volData->getVoxelAt(x+1,y,z).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ - 0.5f), Vector3DFloat(-1.0f, 0.0f, 0.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(-1.0f, 0.0f, 0.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(-1.0f, 0.0f, 0.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(-1.0f, 0.0f, 0.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v1,v2);
						m_meshCurrent->addTriangleCubic(v1,v3,v2);
					}

					int plusYVoxel = m_volData->getVoxelAt(x,y+1,z).getDensity()  >= VoxelType::getThreshold();
					if(currentVoxel > plusYVoxel)
					{
						int material = m_volData->getVoxelAt(x,y,z).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(0.0f, 1.0f, 0.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 1.0f, 0.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(0.0f, 1.0f, 0.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 1.0f, 0.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v1,v2);
						m_meshCurrent->addTriangleCubic(v1,v3,v2);
					}
					if(currentVoxel < plusYVoxel)
					{
						int material = m_volData->getVoxelAt(x,y+1,z).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(0.0f, -1.0f, 0.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, -1.0f, 0.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ - 0.5f), Vector3DFloat(0.0f, -1.0f, 0.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, -1.0f, 0.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v2,v1);
						m_meshCurrent->addTriangleCubic(v1,v2,v3);
					}

					int plusZVoxel = m_volData->getVoxelAt(x,y,z+1).getDensity()  >= VoxelType::getThreshold();
					if(currentVoxel > plusZVoxel)
					{
						int material = m_volData->getVoxelAt(x,y,z).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, 1.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, 1.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, 1.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, 1.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v2,v1);
						m_meshCurrent->addTriangleCubic(v1,v2,v3);
					}
					if(currentVoxel < plusZVoxel)
					{
						int material = m_volData->getVoxelAt(x,y,z+1).getMaterial();

						uint32_t v0 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, -1.0f), material));
						uint32_t v1 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX - 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, -1.0f), material));
						uint32_t v2 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY - 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, -1.0f), material));
						uint32_t v3 = m_meshCurrent->addVertex(PositionMaterialNormal(Vector3DFloat(regX + 0.5f, regY + 0.5f, regZ + 0.5f), Vector3DFloat(0.0f, 0.0f, -1.0f), material));

						m_meshCurrent->addTriangleCubic(v0,v1,v2);
						m_meshCurrent->addTriangleCubic(v1,v3,v2);
					}
				}
			}
		}

		m_meshCurrent->m_Region = m_regSizeInVoxels;

		m_meshCurrent->m_vecLodRecords.clear();
		LodRecord lodRecord;
		lodRecord.beginIndex = 0;
		lodRecord.endIndex = m_meshCurrent->getNoOfIndices();
		m_meshCurrent->m_vecLodRecords.push_back(lodRecord);
	}
}