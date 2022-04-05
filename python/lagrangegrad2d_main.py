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
Nx1=10
Ny2=10

# value of the gradient (/dx,/dy) by the reference triangle rectangle element
# 6 difference possibilities, the first [-1,0] corresponds to the upper-right
# localisation along the center point
derivtrigref = num.array([[-1.,0.],[ 0.,-1.],[1.,-1.],[1.,0.],[0.,1.],[-1.,1.]])


#********** PART 1: IMPORT GEOMETRY FROM inp FILES **********
Ns=[9,2]  # the only parameter that need to be set in hard because it is not in the input file
#that is the number of set ('2' here) and the number of lines to be imported exactly(if not seg. error)

pointsfilename=u'//home//frederic//Documents//github//solidmecanics21//python//untitled3.inp'
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
print("Boundary")
print(str(bcelmts))
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
startbcset=texte.find('*ELSET,ELSET=gamma')
while (not startbcset):
   texte=str(pointsfic.readline())
   startbcset=texte.find('*ELSET,ELSET=gamma')
texte=str(pointsfic.readline())

#bcsetlist={} #dictionnary structure
#bcseti=[]
#iset=0
#t1=[0,0,0,0,0,0,0,0,0,0,0]
#while (iset<=size(Ns,1)):  #number of boundary set
#   while (il<=Ns[iset]):
#      elline=texte.split(',')
#      for i in range(0,size(elline,1)):
#         t1[i]=int(elline[i])
#      texte=str(pointsfic.readline())
#      bcseti.append(t1)
#   bcset{('set'+str(iset))}=bcseti 
#after importing data

#end close the input file
#print("Elements")
#print(str(elmts)) #debug

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
neighbours=[]
for k in range(1,Ne):
   el1=elmts[k] #not elmts[k]
   nn=0
   ng1=[0,0,0]
   i=0
   for l in range(0,Ne):
      if (l==k):
         break;
      #else (elements not identical)
      el2=elmts[l]
      cv=common_vertex(el1,el2)
#      print(str(k)+" "+str(el1)+" : "+str(l)+" "+str(el2)+" : "+str(cv)) #debug
      if (i>=3):
         print("Error of adjacents triangles elements")
         break;
      if ((cv !=-1) and (i<=2)):
         ng1[i]=l
         i=i+1
   neighbours.append(ng1)

#debug
#print(str(neighbours))
#for k in range(0,1):
#    print(str(neighbours[k][0]), str(neighbours([k][1]), str(neighbours([k][2]))))

# Now that mesh is created, the numeric elements of the Laplacian
# problem can be define
Uvar=num.zeros((1,Np)) # the solution vector
rhs=num.zeros((1,Np))

#fills rhs with scalar products over the cells
#for this three functions
def getmatB(j,x1,y2):
   return derivtrigref[j,:]

def getmatJinv(x1,y2):
   return np.array([[1/x1,0.],[0,1/y2]])

def getmatX(x1,y2):
   return np.array([-y2/6.],[x1/6.])


#loop over the rows of rhs (nodes)
#but we will assume there are only regular meshes 
#of trig rectangle and loop over the elements

for i in range(1,Nx1-1):
  for j in range(1,Ny2-1): 
    k=2*(i*Ny2+j)+1
    # if one side has boundary condition, the scalar product
    # must be replaced by the boundary value
    if (neighbours[k][0]==0): #side 0
       rhs[2*k] = const_phie
    elif (neighbours[k][1]==0): #side 2 not previewed
        print("boundary of hypothenusis not previewed")
    elif (neighbours[k][2]==0): #side 1
       rhs[(i+1)*Ny2+j]=const_phie
    else: #normal integration, exact (linear interpol)
      for l in range(0,6): # switch 6 different faces (~hexa circle)
          if (l==0):
             x1=(points[(i+1)*Ny2+j,0]-points[(i)*Ny2+j,0])
             y2=(points[(i+1)*Ny2+j+1,1]-points[(i)*Ny2+j,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          elif (l==1):
             x1=(points[(i+1)*Ny2+j+1,0]-points[(i+1)*Ny2+j,0])
             y2=(points[(i+1)*Ny2+j,1]-points[(i)*Ny2+j,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          elif (l==2):
             x1=(points[(i)*Ny2+j,0]-points[(i)*Ny2+j-1,0])
             y2=(points[(i+1)*Ny2+j,1]-points[(i)*Ny2+j,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          elif (l==3):
             x1=(points[(i)*Ny2+j,0]-points[(i)*Ny2+j-1,0])
             y2=(points[(i)*Ny2+j,1]-points[(i-1)*Ny2+j-1,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          elif (l==4):
             x1=(points[(i)*Ny2+j,0]-points[(i-1)*Ny2+j-1,0])
             y2=(points[(i)*Ny2+j,1]-points[(i-1)*Ny2+j-1,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          elif (l==5):
             x1=(points[(i)*Ny2+j+1,0]-points[(i)*Ny2+j-1,0])
             y2=(points[(i)*Ny2+j,1]-points[(i-1)*Ny2+j-1,1])  
             B=getmatB(l,x1,y2)  # be conscious of the l-parameter here: 
          Jinv=getmatJinv(x1,y2)
          X=getmatX(x1,x2) 
          det=x1*y2/2
          rhs_l = GC*det*num.matprod(num.matprod(B,Jinv), X)
          rhs[i*Ny2+j] = rhs[i*Ny2+j] + rhs_l
     #there should be have a treatment of the points at the sides of the domain which have less than 6 elements or faces

