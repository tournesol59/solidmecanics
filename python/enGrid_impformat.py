# !C:\Python27\python.exe
#project that reads a input file close to an ascii file with comments (symbol --)
#, treats only the uncommented lines which are list of points (3 coordinates of the point in the line)
# which are grouped by 3 lines to form a facet.
# and these facets are written in uml form in an output file

import string
import math
#project_dir_name = u'C:\\Documents and Settings\\HP\\Mes Documents\\CFD_stage\\enGrid\\lidcav_stl\\'
project_dir_name = u'C:\\Users\\Fr\xe9d\xe9ric\\Documents\\CFD_stage\\enGrid\\lidcav_stl\\'
points_fic_name = u'lidcav_geo_sp_ver3.txt'
points_fic_name = project_dir_name + points_fic_name
points_fic = open(points_fic_name, 'r')
#output file:
faces_fic_name = u'lidcav_S3.stl'
faces_fic_name = project_dir_name + faces_fic_name
faces_fic = open(faces_fic_name, 'w')
faces_fic.write('solid ascii\n')
# first lines as comments (header)

texte = points_fic.readline()
comment = string.find(texte, '--')
while comment > 0:
	  # do nothing, next line
	texte=points_fic.readline()
	uncomment = string.find(text, '--')

# func that accepts two 3-tuples of floats and return a tuple as the difference
def calc_vect(p1, p2):
   v=[0.0, 0.0, 0.0]   #init
   for i in range(0,3):
      v[i]=p2[i]-p1[i]
#   print(' 1: {0} {1} {2}\n'.format(repr(p1[0]), repr(p1[1]), repr(p1[2])) )    
#   print(' 2: {0} {1} {2}\n'.format(repr(p2[0]), repr(p2[1]), repr(p2[2])) )
#   print(' v: {0} {1} {2}\n'.format(repr(v[0]), repr(v[1]), repr(v[2]))  )
   return v
# end calc_vect
#func that calculates the cross product of two vectors
def calc_cross(v1, v2):
   v3=[0.0, 0.0, 0.0]   #init
   v3[0]=v1[1]*v2[2]-v1[2]*v2[1]
   v3[1]=v1[2]*v2[0]-v1[0]*v2[2]
   v3[2]=v1[0]*v2[1]-v1[1]*v2[0]
   return v3
#func that calculate a norm
def calc_norm(v):
   n=sqrt( v[0]**2 + v[1]**2 + v[2]**2)
   return n
# func that calculates the normal to a facet of 3 points
def calc_unormal(p1,p2,p3):
   v12=[0.0, 0.0, 0.0]   #init
   v12=calc_vect(p1,p2)
   #print(' v12 {0} {1} {2}\n'.format(repr(v12[0]), repr(v12[1]), repr(v12[2])) )
   v13=[0.0, 0.0, 0.0]   #init
   v13=calc_vect(p1,p3)
#   print(' v13 {0} {1} {2}\n'.format(repr(v13[0]), repr(v13[1]), repr(v13[2])) )
   vx=[0.0, 0.0, 0.0]   #init
   vx=calc_cross(v12,v13)
#  print(' vx: {0} {1} {2}\n'.format(repr(vx[0]), repr(vx[1]), repr(vx[2])) )
   
   if (abs(vx[0]) < 0.001) and (abs(vx[1]) < 0.001):
	   if vx[2] < 0.0:
		   vn=(0,0,-1)
	   else:
		   vn=(0,0,1)
   if (abs(vx[0]) < 0.001) and (abs(vx[2]) < 0.001):
	   if vx[1] < 0.0:
		   vn=(0,-1,0)
	   else: 
		   vn=(0,1,0)
   if (abs(vx[1]) < 0.001) and (abs(vx[2]) < 0.001):
	   if vx[0] < 0.0:
		   vn=(-1,0,0)
	   else:
		   vn=(1,0,0)

   return vn


texte = points_fic.readline()
comment=string.find(texte, '--')
#while (string.rfind(texte, '\n')<>-1):
while (texte<>'\n'):
#   print texte
   if comment<0:
      coord=string.split(texte)
      p1=[0.0, 0.0, 0.0]   #init
      p1[0]=float(coord[0])
      p1[1]=float(coord[1])
      p1[2]=float(coord[2])
   else:
      texte = points_fic.readline()
      comment=string.find(texte, '--')
      continue

   texte = points_fic.readline()
   comment=string.find(texte, '--')
   if comment<0:
      coord=string.split(texte)
      p2=[0.0, 0.0, 0.0]   #init
      p2[0]=float(coord[0])
      p2[1]=float(coord[1])
      p2[2]=float(coord[2])
   else:
      print "Error: three uncommented lines minimum consecutive"
      texte = points_fic.readline()
      comment=string.find(texte, '--')
      continue
   
   texte = points_fic.readline()
   comment=string.find(texte, '--')
   if comment<0:
      coord=string.split(texte)
      p3=[0.0, 0.0, 0.0]   #init
      p3[0]=float(coord[0])
      p3[1]=float(coord[1])
      p3[2]=float(coord[2])
   else:
      print "Error: three uncommented lines minimum consecutive"
      texte = points_fic.readline()
      comment=string.find(texte, '--')
      continue
      
   vn=calc_unormal(p1,p2,p3)
# use of the format function
   print(' facet normal {0} {1} {2}\n'.format(repr(vn[0]), repr(vn[1]), repr(vn[2]) ));
   faces_fic.write(' facet normal {0} {1} {2}\n'.format(repr(vn[0]), repr(vn[1]), repr(vn[2]) ));
   faces_fic.write("  outer loop\n")
   
   faces_fic.write('   vertex {0} {1} {2}\n'.format(repr(p1[0]), repr(p1[1]), repr(p1[2]) ));

   faces_fic.write('   vertex {0} {1} {2}\n'.format(repr(p2[0]), repr(p2[1]), repr(p2[2]) ));

   faces_fic.write('   vertex {0} {1} {2}\n'.format(repr(p3[0]), repr(p3[1]), repr(p3[2]) ));

   faces_fic.write("  endloop\n")
   faces_fic.write(" endfacet\n")

   texte = points_fic.readline()
   comment=string.find(texte, '--')
#end while
faces_fic.write('endsolid\n')
faces_fic.close()
points_fic.close()

ch=raw_input('veuillez saisir une chaine :')
ch
