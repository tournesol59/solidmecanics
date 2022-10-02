#Laplacian computation
import numpy as num
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
import string
import math

#********** PART 0: CONSTANTS ************
EY = 220000.
vu = 0.25
GC = EY/2/(1+vu)
const_phie = 1.0

#constant regular mesh since I did not find a way to completely import the geom
#as user
Nx1=20
Ny2=10

# value of the gradient (/dx,/dy) by the reference triangle rectangle element
# 6 difference possibilities, the first [-1,0] corresponds to the upper-right
# localisation along the center point
derivtrigref = num.array([[-1.,0.],[ 0.,-1.],[1.,-1.],[1.,0.],[0.,1.],[-1.,1.]])
#OTHER WAY: value of gradient in reference triangle element but with form functions
derivtrigrefN_xy = num.array([[-1.,-1.],[1.,0.],[0.,1.]])

def permute_I(i,j,offset): # rectangle inner grid element (inf trig)
    #(i,j) should correspond to the lower left point in the grid. offset as manny points in the boundary
    return num.array([[1 2 3],[offset+i*Ny2+j, offset+i*Ny2+j+1, offset+(i+1)*Ny2+j]])
def permute_II(i,j,offset): # rectangle inner grid element (sup trig)
    return num.array([[1 2 3],[offset+(i+1)*Ny2+j+1, offset+(i+1)*Ny2+j, offset+i*Ny2+j+1]])


#********** PART 1: IMPORT GEOMETRY FROM inp FILES **********
Ns=[2,2]  # the only parameter that need to be set in hard because it is not in the input file
#that is the number of set ('2' here) and the number of lines to be imported exactly(if not seg. error)

pointsfilename=u'//home//frederic//Documents//github//solidmecanics21//python//untitled2.inp'
pointsfic=open(pointsfilename,'r')
startpts=-1
i=0
while (startpts != 0):
   texte=str(pointsfic.readline())
   startpts=texte.find('*NODE');
   i=i+1
#   print(str(i)+" "+str(startpts)+" "+texte) #debug
#print("detect *NODE: "+str(startpts)+" line: "+str(i)) #debug
Np=0
p1=num.zeros((3))
points=num.zeros((1,3))
texte=str(pointsfic.readline())
comment=texte.find('*')
print("point read comment: "+str(comment))
while (comment != 0):   # test with as many lines as triplets coordinates
   Np=Np+1
#   print(str(Np)+" "+texte) #debug
   coord=texte.split(',') 
   coord[3]=(coord[3].split('\n'))[0]
   p1[0]=float(coord[1])
   p1[1]=float(coord[2])
   p1[2]=float(coord[3])   
   points[-1,:]=p1
   texte=str(pointsfic.readline())
   comment=texte.find('*')
   #print("point read comment: "+str(comment)) #debug
print("end count Np :"+str(Np)) #debug
texte=str(pointsfic.readline())# test with as many lines as begin of line elements
startelbc=texte.find('*ELEMENT, type=T3D2')
while (not startelbc):
   texte=str(pointsfic.readline())
   startelbc=texte.find('*ELEMENT, type=T3D2')
#print("point read T3D2: "+str(startelbc)) #debug
Nb=0
t1=[0,0]
bcelmts=[[0,0]]
texte=str(pointsfic.readline())
startelm=texte.find('*ELEMENT, type=CPS3, ELSET=Surface1')
comment=texte.find('*ELEMENT, type=T3D2')
while (startelm != 0):
   texte=str(pointsfic.readline())
   startelm=texte.find('*ELEMENT, type=CPS3, ELSET=Surface1')
   comment=texte.find('*ELEMENT, type=T3D2')  
   while ((comment != 0) and (startelm != 0)):   # test with as many lines as doublet -points reference
      Nb=Nb+1
      coord=texte.split(',')
      coord[2]=(coord[2].split('\n'))[0]
      t1[0]=int(coord[1])
      t1[1]=int(coord[2])
      #t1[2]=int(coord[3]) #out of range
#      print(str(Nb)+" "+str(t1))
      bcelmts.append(t1[:]) # careful when manipulating copy of list!
      texte=str(pointsfic.readline())
      comment=texte.find('*ELEMENT, type=T3D2')
      startelm=texte.find('*ELEMENT, type=CPS3, ELSET=Surface1')
#      print("point read CPS3: "+str(startelm)) #debug
#end while:
#print("Boundary")
#print(str(bcelmts))
while (not startelm):
   texte=str(pointsfic.readline())# test with as many lines as begin of surface elements
   startelm=texte.find('*ELEMENT, type=CPS3, ELSET=Surface1')

Ne=0
t1=[0,0,0];
elmts=[[0,0,0]];
texte=str(pointsfic.readline())
comment=texte.find('*')
eleof=(texte=='')
while ((comment != 0) and (not eleof)):   # test with as many lines as triplet -points reference
   Ne=Ne+1
#   print(str(Ne)+" "+texte) #debug
   coord=texte.split(',')
#   print(str(coord)) #debug
   coord[3]=(coord[3].split('\n'))[0]
   t1[0]=int(coord[1])
   t1[1]=int(coord[2])
   t1[2]=int(coord[3])   
#   print(str(Ne)+" "+str(t1)) #debug
   elmts.append(t1[:]) #careful: see comment to bc above
   texte=str(pointsfic.readline())
   comment=texte.find('*')
   eleof=(texte=='')
#end while
print("end count Ne: "+str(Ne)) #debug
#import boundary set
texte=str(pointsfic.readline())
startbcset=texte.find('*ELSET,ELSET=top')
if (startbcset==-1):
    print("not the right bc card") #debug
    print(texte)
#(try just one side as non zero Dirichlet boundary to begin)
while (startbcset==-1):
   texte=str(pointsfic.readline())
   startbcset=texte.find('*ELSET,ELSET=top')
   print(texte) #debug
texte=str(pointsfic.readline()) #first input of next loop
print(texte) #debug

bcset={} #dictionnary structure
bcsetlog=1
if (bcsetlog==1):
   iset=0
   while ((iset<=len(Ns)) and (startbcset==0)):  #number of boundary set
      bcseti=[]
      il=0
      while (il<Ns[iset]):
         t1=[0,0,0,0,0,0,0,0,0,0,0]
         elline=texte.split(',')
        # elline[9]=(elline[9].split('\n'))[0] #normally not needed
         for i in range(0,10):
           if (elline[i] !=' \n'):
              t1[i]=int(elline[i])
           else:
              break
         texte=str(pointsfic.readline())
         bcseti.append(t1[:])
         il=il+1
      bcset[('set'+str(iset+1))]=bcseti 
      texte=str(pointsfic.readline())
      startbcset=texte.find('*ELSET,ELSET=') 
      iset=iset+1
#after importing data
 
#end close the input file
print("Boundary")
print(str(bcset['set1'])) #debug

pointsfic.close()
#*********** PART 2: DEFINE DATA STRUCT *************
# define if the elements given have common vertex
def common_vertex(el1,el2):
#test vertex 01 of el1
   if (((el2[0]==el1[0]) and (el2[1]==el1[1]))  
      or ((el2[1]==el1[0]) and (el2[0]==el1[1]))
      or ((el2[0]==el1[0]) and (el2[2]==el1[1]))
      or ((el2[2]==el1[0]) and (el2[0]==el1[1]))
      or ((el2[1]==el1[0]) and (el2[2]==el1[1]))
      or ((el2[2]==el1[0]) and (el2[1]==el1[1]))):
      return 0
   elif (((el2[0]==el1[0]) and (el2[1]==el1[2]))  
      or ((el2[1]==el1[0]) and (el2[0]==el1[2]))
      or ((el2[0]==el1[0]) and (el2[2]==el1[2]))
      or ((el2[2]==el1[0]) and (el2[0]==el1[2]))
      or ((el2[1]==el1[0]) and (el2[2]==el1[2]))
      or ((el2[2]==el1[0]) and (el2[1]==el1[2]))):
      return 2
   elif (((el2[0]==el1[1]) and (el2[1]==el1[2]))  
      or ((el2[1]==el1[1]) and (el2[0]==el1[2]))
      or ((el2[0]==el1[1]) and (el2[2]==el1[2]))
      or ((el2[2]==el1[1]) and (el2[0]==el1[2]))
      or ((el2[1]==el1[1]) and (el2[2]==el1[2]))
      or ((el2[2]==el1[1]) and (el2[1]==el1[2]))):
      return 1
   else:
      return -1 # code no common vertex

#create the table of neighbour by parsing the elements matrix-wise of vertices
neighbours=[]
for k in range(1,Ne):
   el1=elmts[k] #not elmts[k]
   ng1=[0,0,0]
   i=0
   l=0
   allfound=0
   while ((l < Ne)and(not allfound)):
      if (l==k):
         l=l+1
         continue;
      else: #(elements not identical)
         el2=elmts[l]
         cv=common_vertex(el1,el2)
         if (k==18):
            print(str(k)+" "+str(el1)+" : "+str(l)+" "+str(el2)+" : "+str(cv)+" ("+str(i)+")") #debug
         if ((cv !=-1) and (i<=2)):
            ng1[i]=l
            i=i+1
      if (i>=3):
         print("Error of adjacents triangles elements")
         allfound=1
      l=l+1
   neighbours.append(ng1[:])
#algo neighbour is ok but useless except we rely on it for boundary conds
#debug
#print(str(neighbours))
#for k in range(0,1):
#    print(str(neighbours[k][0])) #[0]), str(neighbours([k][1]), str(neighbours([k][2]))))

# Now that mesh is created, the numeric elements of the Laplacian
#note that for untitled3.inp the nodes 1-4,5-48 are the boundary nodes
# problem can be define
Uvar=num.zeros((1,Np)) # the solution vector
rhs=num.zeros((1,Np))
Keij=num.zeros((Np,Np))
#fills rhs with scalar products over the cells
#for this three functions
def getmatB(j,x1,y2):
   return derivtrigref[j,:]

def getmatJinv(x1,y2):
   return np.array([[1/x1,0.],[0,1/y2]])

def getmatX(x1,y2):
   return np.array([-y2/6.],[x1/6.]) # comes from the laplacian as the solution
   # of the linear static analysis section indpendent warping problem(and /6 from discrete?)

def getmatKe(x1,y2,Jinv):
    #only the derivative wrt X is present for the moment
   Ke=num.zeros((3,3))
   Ke[0,0]=(derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) 
   Ke[0,1]=(derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) 
   Ke[0,2]=(derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1]*Jinv[1,0]) 
   Ke[1,0]=(derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) 
   Ke[1,1]=(derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) 
   Ke[1,2]=(derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1]*Jinv[1,0]) 
   Ke[2,0]=(derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[0,0]*Jinv[0,0]+derivtrigrefN_xy[0,1]*Jinv[1,0]) 
   Ke[2,1]=(derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1*Jinv[1,0]) *
            (derivtrigrefN_xy[1,0]*Jinv[0,0]+derivtrigrefN_xy[1,1]*Jinv[1,0]) 
   Ke[2,2]=(derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1]*Jinv[1,0]) *
            (derivtrigrefN_xy[2,0]*Jinv[0,0]+derivtrigrefN_xy[2,1]*Jinv[1,0]) 
   return Ke

#loop over the rows of rhs (nodes)
#but we will assume there are only regular meshes 
#of trig rectangle and loop over the elements
x1=points[1,0]-points[0,0] #constant
y2=points[4,1]-points[0,1]
offset_1=2*(Nx+Ny)-1
for (case in range(1,10):  # loop over 9 cases of indexing
  if (case==1):
#left bottom corner
    ni=0
    permute_TI=num.array([[1,2,3],[ni,2*(Nx+Ny)-1,4]])
    permute_TII=num.array([[1,2,3],[2*(Nx+Ny),4,2*(Nx+Ny)-1]])

  elif (case==2):
#left top corner
    ni=4+Nx-1
    permute_TI=num.array([[1,2,3],[ni-1,offset_1+(Nx-2)*(Ny-1),2]])
    permute_TII=num.array([[1,2,3],[ni,2,offset_1+(Nx-2)*(Ny-1)])

  elif (case==3):
#bottom right corner
    ni=2*(Nx)+Ny+1
    permute_TI=num.array([[1,2,3],[ni,3,offset_1+Ny-2]])
    permute_TII=num.array([[1,2,3],[ni-1,offset_1+Ny-2,3]])
   
  elif (case==4):
#right top corner
    ni=offset_1+(Nx-1)*(Ny-1)-1
    permute_TI=num.array([[1,2,3],[ni-1,2+(Nx+Ny),1+(Nx+Ny)]])
    permute_TII=num.array([[1,2,3],[3,1+(Nx+Ny),2+(Nx+Ny)]])

  elif (case==5):
#left side
    for (l in range(1,Nx-1)):
      i=l+3
      permute_TI=num.array([[1,2,3],[i,offset_1+(Ny-1)*(l-1),offset_1+(Ny-1)*l]])
      permute_TII=num.array([[1,2,3],[ offset_1+(Ny-1)*l, l, offset_1+(Ny-1)*(l-1)]])


  elif (case==6):
#right side
     for (l in range(1,Nx-1)):
      i1=offset_1+Ny-1+l
      i2=2*Nx+Ny-l
      permute_TI=num.array([[1,2,3],[i1, i2+1, i1+Ny-1]])
      permute_TII=num.array([[1,2,3],[i2, i1+Ny-1, i2+1]])
  # call sub

  elif (case==7):
#bottom side
     for (l in range(1,Ny-1)):
      i1=offset_1-l+1
      i2=offset_1+l
      permute_TI=num.array([[1,2,3],[i1, i1-1, i2]])
      permute_TII=num.array([[1,2,3],[i2+1, i2, i1-1]])
  # call sub

  elif (case==8):
#top side 
     for (l in range(1,Ny-1)):
      i1=offset_1+(Nx-3)*(Ny-1)+l-1
      i2=Nx+4-1
      permute_TI=num.array([[1,2,3],[i1, i1+1, i2]])
      permute_TII=num.array([[1,2,3],[i2+1, i2, i1+1]])
  # call sub

  elif (case==9):
#interior domain:
     for (i in range(0,Nx-2)):
      for (j in range(0,Ny-2)):
         i1=offset_1+i*(Ny-1)+j
         i2=offset_1+(i+1)*(Ny-1)+j
         permute_TI=num.array([[1,2,3],[i1, i1+1, i2]]) 
         permute_TII=num.array([[1,2,3],[i2+1,i2,i1+1]])

# Jinv=getmatJinv(x1,y2)
# X=getmatX(x1,x2) 
  det=x1*y2/2
  rhs_l = GC*det*num.matprod(num.matprod(B,Jinv), X)
  rhs[i*Ny2+j] = rhs[i*Ny2+j] + rhs_l
     #there should be have a treatment of the points at the sides of the domain which have less than 6 elements or faces

