#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.9.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/sinatz/Desktop/WorkDir/Homogenization/AutomaticUCs/AutomaticCubicUC/MostRecentVersion/Verification')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()
#length of cube
a=5.
#radios of strut
r=1.2
#LengthDiscretization
Lw=.17142857142857142856
#name of file
name='7T1.2R.1714M27Nu'
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Cylinder_1 = geompy.MakeCylinderRH(r, a)
Translation_1 = geompy.MakeTranslation(Cylinder_1, a, 0, 0)
Translation_2 = geompy.MakeTranslation(Cylinder_1, a, a, 0)
Translation_3 = geompy.MakeTranslation(Cylinder_1, 0, a, 0)
Rotation_1 = geompy.MakeRotation(Cylinder_1, OY, 90*math.pi/180.0)
Rotation_2 = geompy.MakeRotation(Cylinder_1, OX, -90*math.pi/180.0)
Translation_4 = geompy.MakeTranslation(Rotation_1, 0, a, 0)
Translation_5 = geompy.MakeTranslation(Rotation_1, 0, a, a)
Translation_6 = geompy.MakeTranslation(Rotation_1, 0, 0, a)
Translation_7 = geompy.MakeTranslation(Rotation_2, a, 0, 0)
Translation_8 = geompy.MakeTranslation(Rotation_2, a, 0, a)
Translation_9 = geompy.MakeTranslation(Rotation_2, 0, 0, a)
Fuse_1 = geompy.MakeFuseList([Cylinder_1, Translation_1, Translation_2, Translation_3, Rotation_1, Rotation_2, Translation_4, Translation_5, Translation_6, Translation_7, Translation_8, Translation_9], True, True)
Box_1 = geompy.MakeBoxDXDYDZ(a, a, a)
Inc = geompy.MakeCommonList([Fuse_1, Box_1], True)
Mat = geompy.MakeCutList(Box_1, [Inc], True)
geompy.TranslateDXDYDZ(Mat, -a/2, -a/2, -a/2)
geompy.TranslateDXDYDZ(Inc, -a/2, -a/2, -a/2)
Composite = geompy.MakePartition([Inc, Mat], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
Inc_1 = geompy.CreateGroup(Composite, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Inc_1, [2])
Mat_1 = geompy.CreateGroup(Composite, geompy.ShapeType["SOLID"])
geompy.UnionIDs(Mat_1, [146])
Group_1 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_1, [158])
Group_2 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_2, [150])
Group_3 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_3, [152])
Group_4 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_4, [102])
Group_5 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_5, [88])
Group_6 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_6, [130])
Group_7 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_7, [148])
Group_8 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_8, [154])
Group_9 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_9, [79])
Group_10 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_10, [52])
Group_11 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_11, [18])
Group_12 = geompy.CreateGroup(Composite, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_12, [156])
[Inc_1, Mat_1, Group_1, Group_2, Group_3, Group_4, Group_5, Group_6, Group_7, Group_8, Group_9, Group_10, Group_11, Group_12] = geompy.GetExistingSubObjects(Composite, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Translation_2, 'Translation_2' )
geompy.addToStudy( Translation_3, 'Translation_3' )
geompy.addToStudy( Rotation_2, 'Rotation_2' )
geompy.addToStudy( Rotation_1, 'Rotation_1' )
geompy.addToStudy( Translation_4, 'Translation_4' )
geompy.addToStudy( Translation_5, 'Translation_5' )
geompy.addToStudy( Translation_6, 'Translation_6' )
geompy.addToStudy( Translation_7, 'Translation_7' )
geompy.addToStudy( Translation_8, 'Translation_8' )
geompy.addToStudy( Translation_9, 'Translation_9' )
geompy.addToStudy( Fuse_1, 'Fuse_1' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Inc, 'Inc' )
geompy.addToStudy( Mat, 'Mat' )
geompy.addToStudy( Composite, 'Composite' )
geompy.addToStudyInFather( Composite, Inc_1, 'Inc' )
geompy.addToStudyInFather( Composite, Mat_1, 'Mat' )
geompy.addToStudyInFather( Composite, Group_1, 'Group_1' )
geompy.addToStudyInFather( Composite, Group_2, 'Group_2' )
geompy.addToStudyInFather( Composite, Group_3, 'Group_3' )
geompy.addToStudyInFather( Composite, Group_4, 'Group_4' )
geompy.addToStudyInFather( Composite, Group_5, 'Group_5' )
geompy.addToStudyInFather( Composite, Group_6, 'Group_6' )
geompy.addToStudyInFather( Composite, Group_7, 'Group_7' )
geompy.addToStudyInFather( Composite, Group_8, 'Group_8' )
geompy.addToStudyInFather( Composite, Group_9, 'Group_9' )
geompy.addToStudyInFather( Composite, Group_10, 'Group_10' )
geompy.addToStudyInFather( Composite, Group_11, 'Group_11' )
geompy.addToStudyInFather( Composite, Group_12, 'Group_12' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Composite,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(Lw,None,1e-07)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
NETGEN_3D = Mesh_1.Tetrahedron()
Inc_2 = Mesh_1.GroupOnGeom(Inc_1,'Inc',SMESH.VOLUME)
Mat_2 = Mesh_1.GroupOnGeom(Mat_1,'Mat',SMESH.VOLUME)
Group_1_1 = Mesh_1.GroupOnGeom(Group_1,'Group_1',SMESH.FACE)
Group_2_1 = Mesh_1.GroupOnGeom(Group_2,'Group_2',SMESH.FACE)
Group_3_1 = Mesh_1.GroupOnGeom(Group_3,'Group_3',SMESH.FACE)
Group_4_1 = Mesh_1.GroupOnGeom(Group_4,'Group_4',SMESH.FACE)
Group_5_1 = Mesh_1.GroupOnGeom(Group_5,'Group_5',SMESH.FACE)
Group_6_1 = Mesh_1.GroupOnGeom(Group_6,'Group_6',SMESH.FACE)
Group_7_1 = Mesh_1.GroupOnGeom(Group_7,'Group_7',SMESH.FACE)
Group_8_1 = Mesh_1.GroupOnGeom(Group_8,'Group_8',SMESH.FACE)
Group_9_1 = Mesh_1.GroupOnGeom(Group_9,'Group_9',SMESH.FACE)
Group_10_1 = Mesh_1.GroupOnGeom(Group_10,'Group_10',SMESH.FACE)
Group_11_1 = Mesh_1.GroupOnGeom(Group_11,'Group_11',SMESH.FACE)
Group_12_1 = Mesh_1.GroupOnGeom(Group_12,'Group_12',SMESH.FACE)
isDone = Mesh_1.Compute()
[ Inc_2, Mat_2, Group_1_1, Group_2_1, Group_3_1, Group_4_1, Group_5_1, Group_6_1, Group_7_1, Group_8_1, Group_9_1, Group_10_1, Group_11_1, Group_12_1 ] = Mesh_1.GetGroups()
Projection_2D = Mesh_1.Projection2D(geom=Group_1)
Source_Face_1 = Projection_2D.SourceFace(Group_7,None,None,None,None,None)
Projection_2D_1 = Mesh_1.Projection2D(geom=Group_2)
Source_Face_2 = Projection_2D_1.SourceFace(Group_8,None,None,None,None,None)
Projection_2D_2 = Mesh_1.Projection2D(geom=Group_3)
Source_Face_3 = Projection_2D_2.SourceFace(Group_12,None,None,None,None,None)
Projection_2D_3 = Mesh_1.Projection2D(geom=Group_4)
Source_Face_4 = Projection_2D_3.SourceFace(Group_9,None,None,None,None,None)
Projection_2D_4 = Mesh_1.Projection2D(geom=Group_5)
Source_Face_5 = Projection_2D_4.SourceFace(Group_10,None,None,None,None,None)
Projection_2D_5 = Mesh_1.Projection2D(geom=Group_6)
Source_Face_6 = Projection_2D_5.SourceFace(Group_11,None,None,None,None,None)
isDone = Mesh_1.Compute()
[ Inc_2, Mat_2, Group_1_1, Group_2_1, Group_3_1, Group_4_1, Group_5_1, Group_6_1, Group_7_1, Group_8_1, Group_9_1, Group_10_1, Group_11_1, Group_12_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'/home/sinatz/Desktop/WorkDir/Homogenization/AutomaticUCs/AutomaticCubicUC/MostRecentVersion/Verification/'+name+'.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
Sub_mesh_1 = Projection_2D.GetSubMesh()
Sub_mesh_2 = Projection_2D_1.GetSubMesh()
Sub_mesh_3 = Projection_2D_2.GetSubMesh()
Sub_mesh_4 = Projection_2D_3.GetSubMesh()
Sub_mesh_5 = Projection_2D_4.GetSubMesh()
Sub_mesh_6 = Projection_2D_5.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Group_8_1, 'Group_8')
smesh.SetName(Group_9_1, 'Group_9')
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Source_Face_1, 'Source Face_1')
smesh.SetName(Source_Face_2, 'Source Face_2')
smesh.SetName(Projection_2D.GetAlgorithm(), 'Projection_2D')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Source_Face_5, 'Source Face_5')
smesh.SetName(Source_Face_6, 'Source Face_6')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Source_Face_3, 'Source Face_3')
smesh.SetName(Group_1_1, 'Group_1')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Source_Face_4, 'Source Face_4')
smesh.SetName(Group_2_1, 'Group_2')
smesh.SetName(Group_3_1, 'Group_3')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Group_4_1, 'Group_4')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Group_5_1, 'Group_5')
smesh.SetName(Group_6_1, 'Group_6')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Group_7_1, 'Group_7')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Group_11_1, 'Group_11')
smesh.SetName(Group_10_1, 'Group_10')
smesh.SetName(Group_12_1, 'Group_12')
smesh.SetName(Mat_2, 'Mat')
smesh.SetName(Inc_2, 'Inc')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
