	 /*
	 --------------------------------------------------------------------

	 Written by Stéphane COMMEND with the assistance
	 of Patricia BOMME, Yves DUBOIS-PELERIN, Dominique
	 EYHERAMENDY, Milan JIRASEK, Andrzej TRUTY and
	 Thomas ZIMMERMANN

	 September 1998

	 Laboratory of Structural and Continuum Mechanics (LSC)
	 Civil Engineering Department (DGC)
	 Swiss Federal Institute of Technology (EPFL)

	 ---------------------------------------------------------------------

	 Modified by : 

	 Vinh Phu NGUYEN,Stephane BORDAS,Cyrille DUNANT, Quang Tri TRUONG
	 and Ravindra DUDDU

	 September 2005

	 to implement the XFEM for cracks, holes,inclusions and biofilms.

	 */

//  MAIN
//  Solves finite element problems using FEM or XFEM 
//
//  The DEBUG option MUST be used (check in file 'debug.def').
//  See also file 'main2.c'.

#include "domain.h"
#include "freestor.h"
#include "compiler.def"

// include for debug ...
#include "flotarry.h" 
#include "node.h" 
#include "element.h" 
#include "geometry_base.h"
#include "geometryentity.h"
#include "vertex.h"
#include "enrichmentitem.h"
#include "enrichmentfunction.h"
#include "cracktip.h"
#include "crackinterior.h"
#include "flotmtrx.h" 
#include <iostream>
#include <vector>
#include <map>
#include "integrationrule.h" 
#include "fei2dtrilin.h"
#include "gausspnt.h"
#include "material.h"
#include "nullmaterial.h"
#include "auxiliaryfield.h"
#include "auxiliaryfield.cpp"
#include "tri_u.h"
#include <fstream> // filestream flux sur les fichiers
 
int main ()
{
 
  Domain* mesh = new Domain() ;

  mesh ->solveFractureMechanicsProblem();
  //mesh ->solveYourself();

  delete mesh ;
  
  /*
  // *********************************************************************
  // DEBUGGING CODE ...
  
  Domain* mesh = new Domain() ;
 
  // CHECK MESH GEOMETRY INTERACTION (THE MOST IMPORTANT PART)

  std::cout<< " Mesh geo interaction ... " << endl;
  mesh->treatMeshGeoInteractionPhase1();
  mesh->treatMeshGeoInteractionPhase2();
  std::cout<< " Finish mesh-geo interaction " << endl ;

  mesh->treatEnrichment() ;
  
  mesh->resolveConflictsInEnrichment();
  mesh->resolveLinearDependencyForEnrichment();

  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 mesh->giveNode(i+1)->printEnrichedNode();
  }

  FILE       *mlbFile;
  char        mlbFileName[15] = "mesh" ; // name of the file
  strcat(mlbFileName,".m");          // extension .m ( Matlab M file)
  mlbFile = fopen(mlbFileName, "w");

  fprintf (mlbFile,"node = [ ");
  double val1 = 0.0, val2 = 0.0;

  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 val1 = mesh->giveNode(i+1)->giveCoordinate(1);
	 val2 = mesh->giveNode(i+1)->giveCoordinate(2);
	 fprintf (mlbFile,"%15.14e %15.14e \n",val1, val2);
  }
  
  fprintf (mlbFile," ];\n ");

  fprintf(mlbFile,"element = [ ");

  Element* elt;
  int numnode;
  Node* node;

  for (size_t e = 0 ; e < mesh->giveNumberOfElements() ; e++)
  {
	 elt = mesh->giveElement(e+1);
	 numnode = elt->giveNumberOfNodes();
	 for (size_t n = 1 ; n <= numnode ; n++)
	 {
		node = elt->giveNode(n);
		fprintf(mlbFile,"%d ",node->giveNumber());   
		if (n==numnode)
		  fprintf(mlbFile,"\n");
	 }
  }

  fprintf(mlbFile,"];\n\n");

  fprintf (mlbFile,"StepEnrichedNodes = [ ");
  double x,y ;
  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 if(mesh->giveNode(i+1)->isStepEnriched())        
	 {
		x = mesh->giveNode(i+1)->giveCoordinate(1);
		y = mesh->giveNode(i+1)->giveCoordinate(2);
		fprintf (mlbFile,"%15.14e %15.14e \n",x,y);
	 }
  }
  fprintf (mlbFile," ];\n ");

  fprintf (mlbFile,"TipEnrichedNodes = [ ");

  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 if(mesh->giveNode(i+1)->isTipEnriched())        
	 {
		x = mesh->giveNode(i+1)->giveCoordinate(1);
		y = mesh->giveNode(i+1)->giveCoordinate(2);
		fprintf (mlbFile,"%15.14e %15.14e \n",x,y);
	 }
  }
  fprintf (mlbFile," ];\n ");

  fprintf (mlbFile,"JunctionEnrichedNodes = [ ");
  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 if(mesh->giveNode(i+1)->isJunctionEnriched())        
	 {
		x = mesh->giveNode(i+1)->giveCoordinate(1);
		y = mesh->giveNode(i+1)->giveCoordinate(2);
		fprintf (mlbFile,"%15.14e %15.14e \n",x,y);
	 }
  }
  fprintf (mlbFile," ];\n ");

  fprintf (mlbFile,"InterfaceEnrichedNodes = [ ");
  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 if(mesh->giveNode(i+1)->isInterfaceEnriched())        
	 {
		x = mesh->giveNode(i+1)->giveCoordinate(1);
		y = mesh->giveNode(i+1)->giveCoordinate(2);
		fprintf (mlbFile,"%15.14e %15.14e \n",x,y);
	 }
  }
  fprintf (mlbFile," ];\n ");

  // compute Gauss points
  for (size_t i = 0 ; i < mesh->giveNumberOfElements() ; i++)
  {
	 elt = mesh->giveElement(i+1);
	 elt->computeGaussPoints();
  }
  //mesh->exportGaussPointsToMatlab(mlbFileName);

  // Export Matlab commands for plotting
  fprintf (mlbFile,"switch size(element,2)\n");
  fprintf (mlbFile,"  case 3 \n");
  fprintf (mlbFile,"  elementType = 'T3';\n");
  fprintf (mlbFile,"  case 6\n");
  fprintf (mlbFile,"  elementType = 'T6';\n");
  fprintf (mlbFile,"  case 4 \n");
  fprintf (mlbFile,"  elementType = 'Q4';\n");
  fprintf (mlbFile,"  case 8 \n");
  fprintf (mlbFile,"  elementType = 'Q8';\n");
  fprintf (mlbFile,"  otherwise \n");
  fprintf (mlbFile,"  error('Grr. Unknown element type');\n");
  fprintf (mlbFile,"end\n\n");

  fprintf (mlbFile,"figure;\n");
  fprintf (mlbFile,"set(gcf,'Color','white');\n");
  fprintf (mlbFile,"plot_mesh(node,element,elementType,'b-');\n");
  
  // plot crack geometry, enriched nodes  

  // crack geometry 
  fprintf (mlbFile,"hold on\n");
  fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
  fprintf (mlbFile,"  crFig = plot(crack(i,1:2:size(crack,2)),crack(i,2:2:size(crack,2)),'r*-');\n");
  fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
  fprintf (mlbFile,"end\n");
 
  // enriched nodes
  fprintf (mlbFile,"h1 = plot(StepEnrichedNodes(:,1),StepEnrichedNodes(:,2),'ro');\n");
  fprintf (mlbFile,"set(h1,'MarkerSize',8);\n");
  fprintf (mlbFile,"h2 = plot(TipEnrichedNodes(:,1),TipEnrichedNodes(:,2),'bs');\n");
  fprintf (mlbFile,"set(h2,'MarkerSize',8);\n");
  fprintf (mlbFile,"h3 = plot(InterfaceEnrichedNodes(:,1),InterfaceEnrichedNodes(:,2),'bs');\n");
  fprintf (mlbFile,"set(h3,'MarkerSize',8);\n");

  // hole or circular inclusion
  fprintf (mlbFile,"r = 1.0 ;\n");   
  fprintf (mlbFile,"xc = 0.0 ;\n");   
  fprintf (mlbFile,"yc = 0.0 ;\n");   
  fprintf (mlbFile,"theta = -pi:0.1:pi;\n");   
  fprintf (mlbFile,"x = xc + r*cos(theta) ;\n");   
  fprintf (mlbFile,"y = yc + r*sin(theta) ;\n"); 
  fprintf (mlbFile,"plot(x,y,'r') ;\n"); 
  fprintf (mlbFile,"axis off ;\n\n");

  fprintf (mlbFile,"figure;\n");
  fprintf (mlbFile,"plot_mesh(node,element,elementType,'b-');\n");
  fprintf (mlbFile,"hold on\n");
  fprintf (mlbFile,"plot(XGP(:,1),XGP(:,2),'rx');\n");

  fclose (mlbFile) ;

  
  /*
  // EXPORT H(X) TO CHECK COMPUTATION OS STEP FUNCTION
  FILE       *mlbFile;
  char        mlbFileName[15] = "stepfunc" ; // name of the file
  strcat(mlbFileName,".m");          // extension .m ( Matlab M file)
  mlbFile = fopen(mlbFileName, "w");

  EnrichmentItem *enr = mesh->giveEnrichmentItem(1);

  std::vector<Element*> *elems = enr->giveElementsInteractWithMe();

  std::vector<EnrichmentFunction*> *step = enr->giveEnrFuncVector();
  (*step)[0]->findActiveEnrichmentItem(enr);

  fprintf (mlbFile,"splitelems = [ ");
  Node *node ;

  for (size_t i = 0 ; i < elems->size() ; i++)
  {
	 size_t numnode = (*elems)[i]->giveNumberOfNodes();
	 for (size_t n = 0 ; n < numnode ; n++)
	 {
		node = (*elems)[i]->giveNode(n+1);
		fprintf(mlbFile,"%d ",node->giveNumber());   
		if(n == numnode-1)
		  fprintf(mlbFile,"\n");
	 }
  }
  
  fprintf (mlbFile," ];\n ");

  fprintf (mlbFile,"stepfunc = [ ");
  double val ;
  
  for (size_t i = 0 ; i < elems->size() ; i++)
  {
	 size_t numnode = (*elems)[i]->giveNumberOfNodes();
	 for (size_t n = 0 ; n < numnode ; n++)
	 {
		node = (*elems)[i]->giveNode(n+1);
		Mu::Point *p = node->makePoint();
		val = (*step)[0]->EvaluateYourSelfAt(p);
		fprintf(mlbFile,"%3.1e ",val);   
		if(n == numnode-1)
		  fprintf(mlbFile,"\n");
	 }
  }
  
  fprintf (mlbFile," ];\n ");
  fprintf (mlbFile,"figure;\n");
  fprintf (mlbFile,"plot_field(node,splitelems,elementType,stepfunc);\n");
  fprintf (mlbFile,"plot_mesh(node,splitelems,elementType,'k-');\n");
  fprintf (mlbFile,"hold on\n");
  fprintf (mlbFile,"for i = 1 : size(crack,1)\n");
  fprintf (mlbFile,"  crFig = plot(crack(1,1:2:size(crack,2)-1),crack(1,2:2:size(crack,2)),'k');\n");
  fprintf (mlbFile,"  set(crFig,'LineWidth',1.5)\n");
  fprintf (mlbFile,"end\n");

  fclose (mlbFile) ;
  */
  /* 
  // EXPORT THE ASYMPTOTIC FUNCTION TO MATLAB TO PLOT

  FILE       *mlbFile;
  char        mlbFileName[10] = "asymp" ; // name of the file
  strcat(mlbFileName,".m");          // extension .m ( Matlab M file)
  mlbFile = fopen(mlbFileName, "w");

  // computed w.r.t the second crack tip
  CrackTip *tip2 = static_cast<CrackTip*>(mesh->giveEnrichmentItem(3));
  std::vector<EnrichmentFunction*> *fns = tip2->giveEnrFuncVector() ; // 4 asymp. functions

  fprintf (mlbFile,"phi2 = [ ");
  double val1 = 0.0, val2 = 0.0, val3 = 0.0, val4 = 0.0;

  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 Node *node = mesh->giveNode(i+1);
	 Mu::Point *p = node->makePoint();
	 (*fns)[0]->findActiveEnrichmentItem(tip2);
	 val1 = (*fns)[0]->EvaluateYourSelfAt(p);
	 (*fns)[1]->findActiveEnrichmentItem(tip2);
	 val2 = (*fns)[1]->EvaluateYourSelfAt(p);
	 (*fns)[2]->findActiveEnrichmentItem(tip2);
	 val3 = (*fns)[2]->EvaluateYourSelfAt(p);
    (*fns)[3]->findActiveEnrichmentItem(tip2);
	 val4 = (*fns)[3]->EvaluateYourSelfAt(p);

	 fprintf (mlbFile,"%8.4e %8.4e %8.4e %8.4e \n",val1, val2,val3,val4);
  }
  
  fprintf (mlbFile," ];\n ");
  /*
  // computed w.r.t the first crack tip
  CrackTip *tip1 = static_cast<CrackTip*>(mesh->giveEnrichmentItem(2));
  std::vector<EnrichmentFunction*> *fns1 = tip1->giveEnrFuncVector() ; // 4 asymp. functions

  fprintf (mlbFile,"phi1 = [ ");
  
  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 Node *node = mesh->giveNode(i+1);
	 Mu::Point *p = node->makePoint();
	 (*fns1)[0]->findActiveEnrichmentItem(tip1);
	 double val1 = (*fns1)[0]->EvaluateYourSelfAt(p);
	 (*fns1)[1]->findActiveEnrichmentItem(tip1);
	 double val2 = (*fns1)[1]->EvaluateYourSelfAt(p);
	 (*fns1)[2]->findActiveEnrichmentItem(tip1);
	 double val3 = (*fns1)[2]->EvaluateYourSelfAt(p);
    (*fns1)[3]->findActiveEnrichmentItem(tip1);
	 double val4 = (*fns1)[3]->EvaluateYourSelfAt(p);

	 fprintf (mlbFile,"%8.4e %8.4e %8.4e %8.4e \n",val1, val2,val3,val4);
  }
  fprintf (mlbFile," ];\n ");

  // write Matlab commands for plotting these functions
  fprintf (mlbFile," subplot(2,2,1);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi2(:,1));;\n ");
  fprintf (mlbFile," title('function B_1');\n ");

  fprintf (mlbFile," subplot(2,2,2);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi2(:,2));;\n ");
  fprintf (mlbFile," title('function B_2');\n ");

  fprintf (mlbFile," subplot(2,2,3);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi2(:,3));;\n ");
  fprintf (mlbFile," title('function B_3');\n ");

  fprintf (mlbFile," subplot(2,2,4);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi2(:,4));;\n ");
  fprintf (mlbFile," title('function B_4');\n ");

  // plot for the first crack tip
  /*
  fprintf (mlbFile," subplot(2,2,1);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi1(:,1));;\n ");
  fprintf (mlbFile," title('function B_1');\n ");

  fprintf (mlbFile," subplot(2,2,2);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi1(:,2));;\n ");
  fprintf (mlbFile," title('function B_2');\n ");

  fprintf (mlbFile," subplot(2,2,3);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi1(:,3));;\n ");
  fprintf (mlbFile," title('function B_3');\n ");

  fprintf (mlbFile," subplot(2,2,4);\n ");
  fprintf (mlbFile," plot_field(node,element,elementType,phi1(:,4));;\n ");
  fprintf (mlbFile," title('function B_4');\n ");

  fclose (mlbFile) ;
  */

  /*
  
  // EXPORT THE STIFFNESS MATRICES TO MATLAB TO CHECK

  FILE       *mlbFile;
  char        mlbFileName[64] = "stiffness" ; // name of the file
  strcat(mlbFileName,".m");          // extension .m ( Matlab M file)
  mlbFile = fopen(mlbFileName, "w");

  //classical elements
  Element *e3 = mesh->giveElement(3);
  FloatMatrix *e3stiffness = e3->giveStiffnessMatrix();
  
  fprintf (mlbFile,"classicalStiffness = [ ");
  
  for (size_t i = 0 ; i < e3stiffness->giveNumberOfRows() ; i++)
  {
	 for (size_t j = 0 ; j < e3stiffness->giveNumberOfColumns() ; j++)
	 {
		double val = e3stiffness->at(i+1,j+1) ;
	   fprintf (mlbFile,"%8.4e ",val);
	 }
	 fprintf (mlbFile,"\n");
  }
  
  fprintf (mlbFile," ];\n ");

  // tip elements
  Element *e163 = mesh->giveElement(163);
  e163->printMyEnrItems();
  FloatMatrix *e163stiffness = e163->giveStiffnessMatrix();

  fprintf (mlbFile,"tip1Stiffness = [ ");
  
  for (size_t i = 0 ; i < e163stiffness->giveNumberOfRows() ; i++)
  {
	 for (size_t j = 0 ; j < e163stiffness->giveNumberOfColumns() ; j++)
	 {
		double val = e163stiffness->at(i+1,j+1) ;
	   fprintf (mlbFile,"%8.4e ",val);
	 }
	 fprintf (mlbFile,"\n");
  }
  
  fprintf (mlbFile," ];\n ");

  Element *e1146 = mesh->giveElement(1146);
  
  FloatMatrix *e1146stiffness = e1146->giveStiffnessMatrix();

  fprintf (mlbFile,"tip2Stiffness = [ ");
  
  for (size_t i = 0 ; i < e1146stiffness->giveNumberOfRows() ; i++)
  {
	 for (size_t j = 0 ; j < e1146stiffness->giveNumberOfColumns() ; j++)
	 {
		double val = e1146stiffness->at(i+1,j+1) ;
	   fprintf (mlbFile,"%8.4e ",val);
	 }
	 fprintf (mlbFile,"\n");
  }
  
  fprintf (mlbFile," ];\n ");

  // split elements
  Element *e484 = mesh->giveElement(495);
  /*
  e484->printMyEnrItems();
  for(size_t i = 0 ; i < e484->giveNumberOfNodes() ; i++)
  {
  Node *noeud = e484->giveNode(i+1) ;
  noeud->printEnrichedNode();
  } */
  /*
  FloatMatrix *e484stiffness = e484->giveStiffnessMatrix();

  fprintf (mlbFile,"splitStiffness = [ ");
  
  for (size_t i = 0 ; i < e484stiffness->giveNumberOfRows() ; i++)
  {
	 for (size_t j = 0 ; j < e484stiffness->giveNumberOfColumns() ; j++)
	 {
		double val = e484stiffness->at(i+1,j+1) ;
	   fprintf (mlbFile,"%8.4e ",val);
	 }
	 fprintf (mlbFile,"\n");
  }
  
  fprintf (mlbFile," ];\n ");

  // plot enriched nodes to check Circle::in(Point*)
  fprintf (mlbFile,"enrichedNodes = [ ");
  double x,y ;
  for (size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
	 if(mesh->giveNode(i+1)->getIsEnriched())         // enriched node
	 {
		x = mesh->giveNode(i+1)->giveCoordinate(1);
		y = mesh->giveNode(i+1)->giveCoordinate(2);
		fprintf (mlbFile,"%15.14e %15.14e \n",x,y);
	 }
  }
  fprintf (mlbFile," ];\n ");

  // check neighbors of a given element. OK :)

  std::vector<Element*> *elems = 
	 static_cast<CrackTip*>(mesh->giveEnrichmentItem(3))->giveElementsInteractWithMe();

  Element *tipElem = (*elems)[0];

  std::list<Element*>* neis = tipElem->giveNeighboringElements();

  fprintf (mlbFile,"neighbors = [ ");
  for(std::list<Element*> :: iterator i = neis->begin() ; i != neis->end() ; i++)
  {
	 Mu::Point *p = (*i)->giveMyCenter();
	 fprintf (mlbFile,"%15.14e %15.14e \n",p->x,p->y);
  }

  fprintf (mlbFile," ];\n ");

  // check Element::in(circle*).Ok
  Mu::Point *tp = new Mu::Point(6,7.5);
  Mu::Circle *c = new Mu::Circle(2,tp);
  if(tipElem->in(c))
	 std::cout << " Ok for in(circle) " << std::endl;

  
  // plot elements belong to the circle ( used to compute SIFs)
  
 
  std::list<Element*> els = tipElem->conflicts(c) ;
  std::cout << " Number of conflicts elements " << els.size() << std::endl;
  fprintf (mlbFile,"intDomain = [ ");
  for(std::list<Element*> :: iterator i = els.begin() ; i != els.end() ; i++)
  {
	 Mu::Point *p = (*i)->giveMyCenter();
	 fprintf (mlbFile,"%15.14e %15.14e \n",p->x,p->y);
  }

  fprintf (mlbFile," ];\n ");
  

  fclose (mlbFile) ;


  /*
  Element *e1738 = mesh->giveElement(1738);
  Mu::Triangle *t = static_cast<Tri_U*>(e1738)->makeTriangle();
  CrackInterior *enr = static_cast<CrackInterior*>(mesh->giveEnrichmentItem(1));
  std::vector<Segment*> *seg = static_cast<PiecewiseLinear*>(enr->giveMyGeo())->giveSegments();
  for(size_t i= 0; i<seg->size();i++)
  {
	 (*seg)[i]->print();
	 if((*seg)[i]->intersects(t)) 
		std::cout << "Ok" ;
  } 

  
  /*

  std::cout << "Finally ..."  << std::endl;
  std::cout << "Set of enriched nodes are :" << std::endl;
  for(size_t i = 0 ; i < mesh->giveNumberOfNodes() ; i++)
  {
  Node *node = mesh->giveNode(i+1) ;
  node->printEnrichedNode();
  }

  */
  /*
  // TEST ON LOCATION ARRAY OF NODE, OF ELEMENT
  Element *e8 = mesh->giveElement(8) ; // 1 node enirched with H(x)
  
  Node *nut;
  for(size_t i = 0 ; i < 3 ; i++)
  {
	 nut = e8->giveNode(i+1);
	 std::cout << "numberof DOF" << nut->giveNumberOfDofs();
	 std::cout << "numberof true DOF" << nut->giveNumberOfTrueDofs();
	 nut->giveLocationArray()->printYourself();
	 nut->giveStandardLocationArray()->printYourself();
	 cout<<endl;
  }
  Node *n = e8->giveNode(1);
  n->giveEnrichedLocationArray()->printYourself();
  IntArray *loc = e8->giveLocationArray();
  loc->printYourself();
  */

  // TEST THE METHOD COMPUTE LOCAL CRACK TIP COORD
  /*
  EnrichmentItem *tip2 = mesh->giveEnrichmentItem(3);

  Mu::Point *p1 = new Mu::Point(7,8.5);
  FloatArray *localCoord = static_cast<CrackTip*>(tip2)->computeLocalCoordOf(p1);
  std::cout << " Local second crack tip coord of P1 : " << std::endl;
  std::cout << " x: " << (*localCoord)[0] << "  y: " << (*localCoord)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord)[1],(*localCoord)[0]) << std::endl;

  Mu::Point *p2 = new Mu::Point(5,8.5);
  FloatArray *localCoord1 = static_cast<CrackTip*>(tip2)->computeLocalCoordOf(p2);
  std::cout << " Local first crack tip coord of P2 : " << std::endl;
  std::cout << " x: " << (*localCoord1)[0] << "  y: " << (*localCoord1)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord1)[1],(*localCoord1)[0]) << std::endl;

  
  Mu::Point *p3 = new Mu::Point(5,6.5);
  FloatArray *localCoord2 = static_cast<CrackTip*>(tip2)->computeLocalCoordOf(p3);
  std::cout << " Local first crack tip coord of P3 : " << std::endl;
  std::cout << " x: " << (*localCoord2)[0] << "  y: " << (*localCoord2)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord2)[1],(*localCoord2)[0]) << std::endl;

  
  Mu::Point *p4 = new Mu::Point(7,6.5);
  FloatArray *localCoord3 = static_cast<CrackTip*>(tip2)->computeLocalCoordOf(p4);
  std::cout << " Local first crack tip coord of P4 : " << std::endl;
  std::cout << " x: " << (*localCoord3)[0] << "  y: " << (*localCoord3)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord3)[1],(*localCoord3)[0]) << std::endl;
  
 
  EnrichmentItem *tip1 = mesh->giveEnrichmentItem(2);

  Mu::Point *p1 = new Mu::Point(3,6.5);
  FloatArray *localCoord1 = static_cast<CrackTip*>(tip1)->computeLocalCoordOf(p1);
  std::cout << " Local first crack tip coord of p1 : " << std::endl;
  std::cout << " x: " << (*localCoord1)[0] << "  y: " << (*localCoord1)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord1)[1],(*localCoord1)[0]) << std::endl;

  Mu::Point *p2 = new Mu::Point(5,6.5);
  FloatArray *localCoord2 = static_cast<CrackTip*>(tip1)->computeLocalCoordOf(p2);
  std::cout << " Local first crack tip coord of p2 : " << std::endl;
  std::cout << " x: " << (*localCoord2)[0] << "  y: " << (*localCoord2)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord2)[1],(*localCoord2)[0]) << std::endl;

  Mu::Point *p3 = new Mu::Point(5,8.5);
  FloatArray *localCoord3 = static_cast<CrackTip*>(tip1)->computeLocalCoordOf(p3);
  std::cout << " Local first crack tip coord of p3 : " << std::endl;
  std::cout << " x: " << (*localCoord3)[0] << "  y: " << (*localCoord3)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord3)[1],(*localCoord3)[0]) << std::endl;
  
  Mu::Point *p4 = new Mu::Point(3,8.5);
  FloatArray *localCoord4 = static_cast<CrackTip*>(tip1)->computeLocalCoordOf(p4);
  std::cout << " Local first crack tip coord of p4 : " << std::endl;
  std::cout << " x: " << (*localCoord4)[0] << "  y: " << (*localCoord4)[1] << std::endl;
  std::cout << " theta: " << atan2((*localCoord4)[1],(*localCoord4)[0]) << std::endl;
  
  */

  /* 
  // TEST OF METHOD GLOBAL2LOCAL() OF FEI2DTRILIN. OKKKKKKK :)
  Element *e59 = mesh->giveElement(59); // a split element
  Mu::Point *p0 = new Mu::Point(e59->giveNode(1)->giveCoordinate(1),e59->giveNode(1)->giveCoordinate(2)) ;
  Mu::Point *p1 = new Mu::Point(e59->giveNode(2)->giveCoordinate(1),e59->giveNode(2)->giveCoordinate(2)) ;
  Mu::Point *p2 = new Mu::Point(e59->giveNode(3)->giveCoordinate(1),e59->giveNode(3)->giveCoordinate(2)) ;

  Mu::Point *localP0 = e59->giveFEInterpolation()->global2Local(mesh,e59->giveNodeArray(),p0);
  localP0->print();

  Mu::Point *localP1 = e59->giveFEInterpolation()->global2Local(mesh,e59->giveNodeArray(),p1);
  localP1->print();

  Mu::Point *localP2 = e59->giveFEInterpolation()->global2Local(mesh,e59->giveNodeArray(),p2);
  localP2->print();

  Mu::Point *midPoint01 = new Mu::Point(0.5*(p0->x+p1->x),0.5*(p0->y+p1->y));
  Mu::Point *localmidPoint01 = e59->giveFEInterpolation()->global2Local(mesh,e59->giveNodeArray(),midPoint01);
  localmidPoint01->print();
  */

  /*
  //TEST ON NUMERICAL INTEGRATION... TAM ON, SE COME BACK TO YOU LATER :)

  Element *e35 = mesh->giveElement(35); // a split element
  e35->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e35->giveNumberOfGaussPoints() << std::endl ;

  GaussPoint **gpArray = e35->giveGaussPointArray();
  for (size_t i = 0 ; i < e35->giveNumberOfGaussPoints() ; i++) {
  GaussPoint *gp = gpArray[i];
  Mu::Point *p = gp->giveCoordinates();
  std::cout << " Coordinate of Gauss points :" << std::endl;
  p->print();
  }

  Element *e186 = mesh->giveElement(186); // a split element
  //e37->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e186->giveNumberOfGaussPoints() << std::endl ;
  
  /*
  Element *e555 = mesh->giveElement(555); // a split element
  //e555->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e555->giveNumberOfGaussPoints() << std::endl ;

  Element *e467 = mesh->giveElement(467); // a split element
  //e467->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e467->giveNumberOfGaussPoints() << std::endl ;

  Element *e13 = mesh->giveElement(13); 
  e13->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e13->giveNumberOfGaussPoints() << std::endl ;

  GaussPoint **gaussPointArray = e13->giveGaussPointArray();
  for (size_t i = 0 ; i < e13->giveNumberOfGaussPoints() ; i++) {
  GaussPoint *gp = gaussPointArray[i];
  Mu::Point *p = gp->giveCoordinates();
  std::cout << " Coordinate of Gauss points :" << std::endl;
  p->print();
  }


  // Tip elements

  Element *e361 = mesh->giveElement(361); // a tip element
  e361->computeGaussPoints();
  std::cout << "Number of Gauss point : " ;
  std::cout << e361->giveNumberOfGaussPoints() << std::endl ;

  */

  /*
  // test the stiffnes matrix of split element
 
  FloatMatrix *ma35 = e35->computeStiffnessMatrix();

  for(size_t i = 0 ; i < e35->giveNumberOfNodes() ; i++)
  {
  Node *noeud = e35->giveNode(i+1) ;
  noeud->printEnrichedNode();
  }
  // test the size of the enriched matrix. OK :)
  std::cout << ma35->giveNumberOfColumns() << std::endl ;


  /*
  // TEST OF THE AUXILIARY FIELD CALCULATION 

  AuxiliaryField<Material,NullMaterial,Mode_i,PlaneStrain> * auxFieldHomo 
  = new AuxiliaryField<Material,NullMaterial,Mode_i,PlaneStrain>();

  CrackTip *tip = static_cast<CrackTip*>(mesh->giveEnrichmentItem(3)); // a tip
  Mu::Point *p = new Mu::Point(3,4) ;

  FloatMatrix *dU2dx = auxFieldHomo->ComputedU2dx(tip,p);  

  Element *ele = mesh->giveElement(3);
  FloatMatrix *Dmatrix = ele->giveConstitutiveMatrix();
  Dmatrix->printYourself();

  FloatArray *strain = auxFieldHomo->ComputeStrainVector(dU2dx); 
  FloatArray *stress = auxFieldHomo->ComputeStressVector(Dmatrix,strain);  
  stress->printYourself();

  std::cout << std::endl ;

  tip->printYourSelf(); */

  // test of computation of enrichment functions
  /*

  EnrichmentItem *tip2 = mesh->giveEnrichmentItem(3);
  tip2->printYourSelf();
  vector<EnrichmentFunction*> *myEnr = tip2->giveEnrFuncVector();
  Mu::Point *p = new Mu::Point(7.0,8.5);
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function :" << std::endl;
  (*myEnr)[i]->findActiveEnrichmentItem(tip2);
  std::cout << (*myEnr)[i]->EvaluateYourSelfAt(p)<< std::endl;
  }
  FloatArray *ret ;
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function derivative :" << std::endl;
  ret = (*myEnr)[i]->EvaluateYourGradAt(p);
  ret->printYourself();
  }
  */

  /*
  EnrichmentItem *tip1 = mesh->giveEnrichmentItem(2);
  tip1->printYourSelf();
  vector<EnrichmentFunction*> *myEnr = tip1->giveEnrFuncVector();
  Mu::Point *p = new Mu::Point(3.0,6.5);
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function :" << std::endl;
  (*myEnr)[i]->findActiveEnrichmentItem(tip1);
  std::cout << (*myEnr)[i]->EvaluateYourSelfAt(p)<< std::endl;
  }
  FloatArray *ret ;
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function derivative :" << std::endl;
  ret = (*myEnr)[i]->EvaluateYourGradAt(p);
  ret->printYourself();
  }

  EnrichmentItem *interior = mesh->giveEnrichmentItem(1);
  interior->printYourSelf();
  vector<EnrichmentFunction*> *myEnr = interior->giveEnrFuncVector();
  Mu::Point *p = new Mu::Point(5.0,5.5);
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function :" << std::endl;
  (*myEnr)[i]->findActiveEnrichmentItem(interior);
  std::cout << (*myEnr)[i]->EvaluateYourSelfAt(p)<< std::endl;
  }

  FloatArray *ret ;
  for(size_t i = 0 ; i < myEnr->size() ; i++)
  {
  (*myEnr)[i]->printMyEnrichmentItem();std::cout << std::endl;
  std::cout << " Value of the " << i+1 << " th enrichment function derivative :" << std::endl;
  ret = (*myEnr)[i]->EvaluateYourGradAt(p);
  ret->printYourself();
  }
  */
  /*
  Element *e = mesh->giveElement(1);
  Mu::Point *p = new Mu::Point(1,0);
  Mu::Point *glo
  = e->giveFEInterpolation()->local2Global(mesh,e->giveNodeArray(),p);
  glo->print();
  */
  
  /*
  // TEST ON ELEMENT::CONFLICTS(CIRCLE*)

  Vertex *tip = static_cast<Vertex*>(mesh->giveGeoEntity(2));
  double x = tip->giveCoordinate(1);
  double y = tip->giveCoordinate(2);
  double r = 4.0 ;


  Mu::Point *center = new Mu::Point(x,y);


  Element *e1146 = mesh->giveElement(1146);
  double A = e1146->area();
  std::cout<< " Area of element 1146 using Element::area() "<< A << std::endl;
  Mu::Triangle *t = static_cast<Tri_U*>(e1146)->makeTriangle();
  std::cout<< " Area of element 1146 using Triangle::area() "<< t->area() << std::endl;

  Mu::Circle *ball  = new Mu::Circle(r*A,center);
  std::cout << " Element 1146 belongs to the ball ? " ;
  if(e1146->in(ball))
  std::cout << " YES :) " << std::endl;
  else
  std::cout << " No, WHY???" << std::endl ;

  Element *e1 = mesh->giveElement(1);
  std::cout << " Element 1 belongs to the ball ? " ;
  if(e1->in(ball))
  std::cout << " YES :) "<< std::endl;
  else
  std::cout << " No, WHY??? "<< std::endl;


  Element *e101 = mesh->giveElement(101);
  std::cout << " Element 101 belongs to the ball ? " ;
  if(e101->in(ball))
  std::cout << " YES :) "<< std::endl;
  else
  std::cout << " No, WHY??? "<< std::endl;

  Element *e667 = mesh->giveElement(667);
  std::cout << " Element 667 belongs to the ball ? " ;
  if(e667->in(ball))
  std::cout << " YES :) "<< std::endl;
  else
  std::cout << " No, WHY??? "<< std::endl;

  std::list<Element*> dom = e1146->conflicts(ball);
  for( std::list<Element*>::iterator it = dom.begin() ; it != dom.end() ; it++)
  std::cout << (*it)->giveNumber() << std::endl ;
  
  
  for(size_t i = 0 ; i < mesh->giveNumberOfElements() ; i++)
	 mesh->giveElement(i+1)->clearChecked();


  // tip 2
  Vertex *tip2 = static_cast<Vertex*>(mesh->giveGeoEntity(3));


  double x2 = tip2->giveCoordinate(1);
  double y2 = tip2->giveCoordinate(2);
  double r2 = 4.0 ;


  Mu::Point *center2 = new Mu::Point(x2,y2);


  Element *e163 = mesh->giveElement(163);
  e163->printMyEnrItems();
  double A2 = e163->area();
  std::cout<< " Area of element 1146 using Element::area() "<< A2 << std::endl;
  Mu::Triangle *t2 = static_cast<Tri_U*>(e163)->makeTriangle();
  std::cout<< " Area of element 163 using Triangle::area() "<< t2->area() << std::endl;

  Mu::Circle *ball2  = new Mu::Circle(r2*A2,center2);

  std::cout << " Element 163 belongs to the ball ? " ;
  if(e163->in(ball2))
	 std::cout << " YES :) " << std::endl;
  else
	 std::cout << " No, WHY???" << std::endl ;

  std::list<Element*> dom2 = e163->conflicts(ball2);
  for( std::list<Element*>::iterator it = dom2.begin() ; it != dom.end() ; it++)
	 std::cout << (*it)->giveNumber() << std::endl ;
  */
  /*
  // TEST OF COMPUTATION OF BRANCH FUNCTIONS
  EnrichmentItem *crackTip2 = mesh->giveEnrichmentItem(3);
  std::vector<EnrichmentFunction*>* enrFns = crackTip2->giveEnrFuncVector();
  Mu::Point *globalCoordGP = new Mu::Point(7,8.5);
  for(size_t i = 0 ; i < enrFns->size() ; i++)
  {
	 EnrichmentFunction* enrFn = (*enrFns)[i];
	 enrFn->findActiveEnrichmentItem(crackTip2);
	 double phiGP = enrFn->EvaluateYourSelfAt(globalCoordGP);
	 std::cout << " Function Phi " << i+1 << "th : " << phiGP << std::endl ;
	 FloatArray *gradPhiGP = enrFn->EvaluateYourGradAt(globalCoordGP);

	 double dPhidXGP = (*gradPhiGP)[0] ;
	 std::cout << " dPhidx of " << i+1 << "th function : " << dPhidXGP << std::endl ;
	 double dPhidYGP = (*gradPhiGP)[1] ; 		
	 std::cout << " dPhidy of " << i+1 << "th function : " << dPhidYGP << std::endl ;
  }
  */
  /*
  // TEST ON METHOD INTERSECTS 
  Mu::Segment *s1 = new Mu::Segment(Mu::Point(1,0.5),Mu::Point(3,0.5));
  Mu::Segment *s2 = new Mu::Segment(Mu::Point(-1.,0.5),Mu::Point(2.,0.5));
  Mu::Triangle *tri = new Mu::Triangle(&Mu::Point(1.,0.),&Mu::Point(2.,0.),&Mu::Point(1.,1.));
  
  if(s1->intersects(tri))
	 std::cout << "OK" << std::endl ;
  else
    std::cout << " not OK" ;

  std::vector<Mu::Point> inte = s1->intersection(tri);
  std::cout << "there are " << inte.size() << " intersection points " << std::endl ;

  if(s2->intersects(tri))
    std::cout << "OK" << std::endl ;
  else
    std::cout << " not OK" ;

  std::vector<Mu::Point> gd = s2->intersection(tri);
  std::cout << "there are " << gd.size() << " intersection points " << std::endl ;

  Mu::Segment *s3 = new Mu::Segment(Mu::Point(0,0),Mu::Point(1,1));
  Mu::Point *pp = new Mu::Point(0.5,0.6) ; 
  if (s3->on(pp))
	 std::cout << "pp on this segment . Good :) " << std::endl ;
  else
    std::cout << " pp not on this segment " << std::endl ;

  Mu::Segment *s4 = new Mu::Segment(Mu::Point(0,0),Mu::Point(0,1));
  Mu::Point *ppp = new Mu::Point(0,3) ; 
  if (s4->on(ppp))
	 std::cout << "ppp on this segment . Good :) " << std::endl ;
  else
    std::cout << " ppp not on this segment " << std::endl ;

  // case if segment touches one vertex of the triangle !!!
  Mu::Segment *s5 = new Mu::Segment(Mu::Point(1,0.9),Mu::Point(2,0.9));

  if(s5->intersects(tri))
	 std::cout << "OK" << std::endl ;
  else
    std::cout << " not OK" << std::endl ;

	 Mu::Segment *s1 = new Mu::Segment(Mu::Point(0.,0.),Mu::Point(1.,1.));
  Mu::Segment *s2 = new Mu::Segment(Mu::Point(1.,0.),Mu::Point(0.,1.));
  Mu::Triangle *tri = new Mu::Triangle(&Mu::Point(0.,0.),&Mu::Point(1.,0.),&Mu::Point(0.,1.));
  std::cout<< tri->area() << std::endl;
  s1->print();
  if(s1->intersects(tri))
	 std::cout << "OK" ;
  else
	 std::cout << " not OK" ;

  if(s1->intersects(s2))
	 std::cout << "OK" ;
  else
	 std::cout << " not OK" ;

  Mu::Point p = s1->intersection(s2) ;
  p.print();

  Element *e348 = mesh->giveElement(109);
  Mu::Triangle *t = static_cast<Tri_U*>(e348)->makeTriangle();
  GeometryEntity *geo = mesh->giveGeoEntity(1);
  
  
  if(geo->intersects(t))
	 std::cout << "OK" ;
  else
	 std::cout << " not OK" ;
  */
  /*
  Mu::Point p1(-1,-1);
  Mu::Point p2(1,1);
  Mu::Segment s(p1,p2); 
 
  Mu::Point *pp3=new Mu::Point(0,0);
  Mu::Point *pp4=new Mu::Point(1,0);
  Mu::Point *pp5=new Mu::Point(0,1);

  Mu::Triangle *t = new Mu::Triangle(pp3,pp4,pp5);

  std::vector<Mu::Point> intersects = s.intersection(t);

  if(intersects.size() == 1)
  {
	 std::cout << " Need to check corners of triangle " << std::endl ;
	 if(s.on(t->getBoundingPoint(0)))
		intersects.push_back(*t->getBoundingPoint(0));
	 if(s.on(t->getBoundingPoint(1)))
		intersects.push_back(*t->getBoundingPoint(1));
	 if(s.on(t->getBoundingPoint(2)))
		intersects.push_back(*t->getBoundingPoint(2));
  }
  if(intersects.size() == 2)
    std::cout << " em dep qua !!!  " << std::endl ;

  for(size_t i = 0 ; i < intersects.size() ; i++)
	 intersects[i].print();

  // test method Segment::on(Point*) of Cyrille
  Mu::Segment ss(Mu::Point(-3,0),Mu::Point(-1,0));
  ss.vector()->print();
  Mu::Point *testpoint = new Mu::Point(-2,0);

  if(ss.on(testpoint))
    std::cout << " Greatttttt !!!  " << std::endl ;
  else
    std::cout << " Whyyyyyyy  ???  " << std::endl ;

  if(!isAligned(testpoint, ss.first(),ss.second()))
	 std::cout << " not aligned,  stop " << std::endl ;

  Mu::Segment sss(Mu::Point(3,0),Mu::Point(1,0));
  
  Mu::Point *tp = new Mu::Point(2,0);

  if(sss.on(tp))
    std::cout << " Greatttttt !!!  " << std::endl ;
  */
/*
// Test the allowable tolerance of Segment::on(Point*)

Mu::Point p1(0,1.0001);
Mu::Point p2(2,1.0001);
double y = 1.00011 ;
y = floor(y*1e10+0.5)/1e10;
printf("%3f12",y);
Mu::Point *testPoint = new Mu::Point(1,y);

Mu::Segment s(p1,p2);
if(s.on(testPoint))
  std::cout << "OK" ;
else
  std::cout << "not OK" ;
*/
/*
  // check geometry facilities for multi-segments cracks
  std::cout<< " Testing on geometry facilities of PiecewiseLinear ... " << endl ;

  std::cout<< "  1. Testing on finding of closest segment ... " << endl ;
  EnrichmentItem *crack = mesh->giveEnrichmentItem(1);

  Mu::Point *p1 = new Mu::Point(0.1,1.5); // first segment
  Mu::Point *p2 = new Mu::Point(0.3,1.5); // second segment
  Mu::Point *p3 = new Mu::Point(0.3,0.8); // second segment
  Mu::Point *p4 = new Mu::Point(0.3,0.8); // second segment

  Mu::Segment s1 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->FindSegmentClosestTo(p1);
  Mu::Segment s2 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->FindSegmentClosestTo(p2);
  s1.print();
  s2.print();

  std::cout<< "  2. Testing on signed distance computation ... " << endl ;
  double d1 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->computeSignedDistanceOfPoint(p1);
  double d2 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->computeSignedDistanceOfPoint(p2);
  std::cout << d1 << endl;
  std::cout << d2 << endl;

  std::cout<< "  3. Testing on position of point w.r.t crack ... " << endl ;
  int pos1 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->givePositionComparedTo(p1);
  int pos2 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->givePositionComparedTo(p2);
  int pos3 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->givePositionComparedTo(p3);
  int pos4 = 
	 (static_cast<PiecewiseLinear*>(crack->giveMyGeo()))->givePositionComparedTo(p4);
  std::cout << pos1 << endl;
  std::cout << pos2 << endl;
  std::cout << pos3 << endl;
  std::cout << pos4 << endl;
  */
//char r;
//std :: cin >> r;

return 0;

}

