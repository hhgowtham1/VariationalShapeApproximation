#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <ostream>
#include <stdlib.h> // in C, for rand()
#include <time.h>    //for time()
#include <queue> //priority queue
#include <igl/readOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/facet_adjacency_matrix.h>
#include <vector>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/convex_hull_3.h>
//#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_3.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/unique_rows.h>

#include "HalfedgeBuilder.cpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
//typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt> Delaunay;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3   Point;
//typedef CGAL::Point_3<CGAL::Epick> point;
using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;


MatrixXd V;
MatrixXi F;


MatrixXd N_faces;   //computed calling pre-defined functions of LibiGL
MatrixXd N_vertices; //computed calling pre-defined functions of LibiGL


MatrixXd lib_N_vertices;  //computed using face-vertex structure of LibiGL
MatrixXi lib_Deg_vertices;//computed using face-vertex structure of LibiGL

MatrixXd he_N_vertices; //computed using the HalfEdge data structure
MatrixXi regions;// = MatrixXi::Zero(F.rows(),3);
vector<vector<int>> neighbTri;//= MatrixXi::Zero(F.rows(),3);
MatrixXi regTri;
int k;//=;160;
double threshold;// = 0.00001;
double PI= acos(-1);
int l2orl21;// = 1;
vector<vector<int>> clusters; // no of proxied regions or clusters
vector<int> initialTriangles;
vector<double> errs;
HalfedgeDS* he;
vector<vector<int>> polygon_vertices;



bool vector_contains(vector<int> veclist, int a){
    return count(veclist.begin(), veclist.end(), a);
}
bool vector_contains(vector<double> veclist, double a){
    return count(veclist.begin(), veclist.end(), a);
}


RowVector3d barycenter(Vector3d v1, Vector3d v2, Vector3d v3){
    return (v1+v2+v3)/(3.0);
}

RowVector3d barycenter(MatrixXd V){
    return (V.row(0)+V.row(1)+V.row(2))/(3.0);
}

double areaOfTriangle(Vector3d v1, Vector3d v2, Vector3d v3){
    
    return ((v2-v1).cross(v3-v1)).norm() * 0.5 ;
}


double areaOfTriangle(MatrixXd V){
    RowVector3d v1; RowVector3d v2; RowVector3d v3;
    v1=V.row(0);
    v2=V.row(1);
    v3=V.row(2);
        
//    double area = ((V.row(1)-V.row(0)).cross(V.row(2)-V.row(0))).norm();
//    return area*0.5 ;
    return ((v2-v1).cross(v3-v1)).norm() * 0.5 ;
}

Vector3d triangleNormal(Vector3d v1, Vector3d v2, Vector3d v3){
    return (  (v2-v1).cross(v3-v1) / (  (v2-v1).cross(v3-v1).norm() )  ).normalized();
}

Vector3d triangleNormal(MatrixXd t){
    Vector3d v1 = t.row(0);
    Vector3d v2 = t.row(1);
    Vector3d v3 = t.row(2);
    return (  (v2-v1).cross(v3-v1) / (  (v2-v1).cross(v3-v1).norm() )  ).normalized();
}

Vector3d computeXi(vector<MatrixXd> t){
    RowVector3d Xi=Vector3d(0.0,0.0,0.0);
    double denom=0.0;
    for (int i=0; i<t.size(); i++) {
        Xi += barycenter(t.at(i)) * areaOfTriangle(t.at(i));
        denom += areaOfTriangle(t.at(i));
    }
    return Xi/denom;
}
Vector3d computeXi(vector<int> tri )//for region
{
    vector<MatrixXd> tr ;
    MatrixXd m=MatrixXd::Zero(3,3);
    
    for(int i=0; i<tri.size(); i++){
        m.row(0)=V.row(F(tri.at(i),0));
        m.row(1)=V.row(F(tri.at(i),1));
        m.row(2)=V.row(F(tri.at(i),2));
//        cout<<"\n"<<m<<"\n";
        
        tr.push_back(m);
    }
//    cout<<"going here";
   return computeXi(tr);
}
Vector3d computeXi(MatrixXd t){
    
    vector<MatrixXd> tr;
    tr.push_back(t);
    
    return computeXi(tr);
}
Vector3d computeXi(int t){
    MatrixXd tr=MatrixXd::Zero(3,3);
    for( int i=0;i<3;i++)
    tr.row(i)= V.row(F(t,i));
    return computeXi(tr);
}
Vector3d computel2Ni(vector<MatrixXd> t){
    MatrixXd Mi= MatrixXd::Zero(3,3);
    MatrixXd Cov = MatrixXd::Zero(3,3);
    Vector3d Ni;//= Vector3d::Zero(1,3);
    Vector3d Xi = computeXi(t);
    MatrixXd s = MatrixXd::Zero(3,3);
    s(0,0)= 10.0;   s(0,1)= 7.0; s(0,2) = 0.0;
    s(1,0)= 7.0;    s(1,1)= 10.0; s(1,2) = 0.0;
    s(2,0) = 0.0; s(2,1) = 0.0; s(2,2) = 0.0;
    double secTerm = 0.0;
    for (int i=0; i<t.size(); i++) {
        double area = areaOfTriangle(t.at(i));
        secTerm += area;
        Mi.row(0)= (t.at(i)).row(1)-t.at(i).row(0);
        Mi.row(1)= t.at(i).row(2)-t.at(i).row(0);
        //type cast (double)?
        Cov+= ( ( (2.0 * area / 72.0) * Mi * s * Mi.transpose() ) + area * barycenter(t.at(i)).transpose() * barycenter(t.at(i))    );// -  ( area * computeXi(t.at(i)) * computeXi(t.at(i)).transpose()    ) ;
    }
    Cov = Cov - secTerm* Xi * Xi.transpose();
    EigenSolver<MatrixXd> ESolver(Cov);
    VectorXd Evalues= ESolver.eigenvalues().real();
    MatrixXd Evectors= ESolver.eigenvectors().real();
    
    double minVal;
    MatrixXd::Index minRow, minCol;
    minVal = Evalues.minCoeff(&minRow, &minCol);
//        std::cout<<minindex<<std::endl;
//    Ni= Evectors.col(minindex).transpose();
    Ni = Evectors.col(minRow);
    Ni=Ni.normalized();
    return Ni;
}
Vector3d computel2Ni(MatrixXd t){
    
    return triangleNormal(t);
}

Vector3d computel2Ni(int t){
    MatrixXd tr=MatrixXd::Zero(3,3);
    for( int i=0;i<3;i++)
    tr.row(i)= V.row(F(t,i));
    return triangleNormal(tr);
}



Vector3d computel2Ni(vector<int> t)
{
    vector<MatrixXd> tr;
//    cout<<"size: "<<t.size()<<"\n";
    for(int i=0; i<t.size(); i++)
    {
        MatrixXd a=MatrixXd::Zero(3,3);
        for(int j=0; j<3; j++)
        {
            a.row(j)= V.row(F(t.at(i),j));
        }
        tr.push_back(a);
    }
//    cout<<"goingheralso";
        return computel2Ni(tr);
}
double orthoDistance(Vector3d Xi, Vector3d Ni, Vector3d v){
    Ni=Ni.normalized();
    return fabs((Ni.dot(v-Xi))/Ni.norm());
}



double l2error(Vector3d Xi, Vector3d Ni, MatrixXd T){
    
    double d1,d2,d3;
    d1= orthoDistance(Xi,Ni,T.row(0));
    d2= orthoDistance(Xi,Ni,T.row(1));
    d3= orthoDistance(Xi,Ni,T.row(2));

    return (1.0/6.0) * (d1*d1 + d2*d2 + d3*d3 + d1*d2 + d1*d3 + d2*d3)* areaOfTriangle(T);
}
double l2error(Vector3d Xi, Vector3d Ni, int t){
    MatrixXd tr=MatrixXd::Zero(3,3);
    for( int i=0;i<3;i++)
    tr.row(i)= V.row(F(t,i));
    return l2error(Xi,Ni,tr);
}

double l2error(vector<int> tri, int f){
    MatrixXd t = MatrixXd::Zero(3,3);
    t.row(0)= V.row(F(f,0));
    t.row(1)= V.row(F(f,1));
    t.row(2)= V.row(F(f,2));
    
    vector<MatrixXd> tr ;
    MatrixXd m=MatrixXd::Zero(3,3);
    for(int i=0; i<tri.size(); i++){
        m.row(0)=V.row(F(tri.at(i),0));
        m.row(1)=V.row(F(tri.at(i),1));
        m.row(2)=V.row(F(tri.at(i),2));
        tr.push_back(m);
    }
    return l2error(computeXi(tr) , computel2Ni(tr) , t);
}

double l2error(int t, int f ){
    vector<int> tri;
    tri.push_back(t);
    return l2error(tri, f);
    
}
double l2error(int f){
    
    MatrixXd t = MatrixXd::Zero(3,3);
    t.row(0)= V.row(F(f,0));
    t.row(1)= V.row(F(f,1));
    t.row(2)= V.row(F(f,2));
    
    return l2error(computeXi(t), computel2Ni(t), t);
}

Vector3d computel21Ni(vector<MatrixXd> tr)
{
    Vector3d Ni;
    for(int i=0; i<tr.size(); i++)
    {
        Ni += areaOfTriangle(tr[i]) * triangleNormal(tr[i]) ;
    }
    Ni= Ni.normalized();
//    cout<<"N:= "<<Ni.transpose()<<endl;
    return Ni;
}

Vector3d computel21Ni(MatrixXd tr)
{
    vector<MatrixXd> t;
    t.push_back(tr);
    return computel21Ni(t);
}

Vector3d computel21Ni(vector<int> t)
{
    vector<MatrixXd> tr;
//    cout<<"size: "<<t.size()<<"\n";
    for(int i=0; i<t.size(); i++)
    {
        MatrixXd a=MatrixXd::Zero(3,3);
        for(int j=0; j<3; j++)
        {
            a.row(j)= V.row(F(t.at(i),j));
        }
        tr.push_back(a);
    }
//    cout<<"goingheralso";
        return computel21Ni(tr);
}

double l21error(Vector3d Ni, MatrixXd T){
    Ni=Ni.normalized();
    return ((triangleNormal(T) - Ni).squaredNorm()) * areaOfTriangle(T);
}
double l21error(Vector3d Ni, int t){
//    cout<<"N: "<<Ni.transpose()<<endl;
    Ni=Ni.normalized();
    MatrixXd tr=MatrixXd::Zero(3,3);
    for( int i=0;i<3;i++)
    tr.row(i)= V.row(F(t,i));
    return l21error(Ni,tr);
}
double l21error(int f){
    
    MatrixXd t = MatrixXd::Zero(3,3);
    t.row(0)= V.row(F(f,0));
    t.row(1)= V.row(F(f,1));
    t.row(2)= V.row(F(f,2));
    
    return l21error( computel21Ni(t), t);
}
double l21error(vector<int> tri, int f){
    MatrixXd t = MatrixXd::Zero(3,3);
    t.row(0)= V.row(F(f,0));
    t.row(1)= V.row(F(f,1));
    t.row(2)= V.row(F(f,2));
    
    vector<MatrixXd> tr ;
    MatrixXd m=MatrixXd::Zero(3,3);
    for(int i=0; i<tri.size(); i++){
        m.row(0)=V.row(F(tri.at(i),0));
        m.row(1)=V.row(F(tri.at(i),1));
        m.row(2)=V.row(F(tri.at(i),2));
        tr.push_back(m);
    }
    Vector3d N= computel21Ni(tr);
//    cout<<"N: "<<N.transpose()<<endl;
    return l21error(N , t);
}

struct for_region {
    int face;
    int cluster;
    double error;
    for_region(int face, int cluster, double error) : face(face), cluster(cluster), error(error){}
};
struct CompareError {
    bool operator()(for_region const& t1, for_region const& t2)
    {
        return t1.error < t2.error;
    }
};

double disrtortionError(vector<vector<int>> clusters)
{

    MatrixXi R(F.rows(),1);
    for(int i=0; i<F.rows(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(vector_contains(clusters.at(j),i))
                R(i,0) = j;
    }

    double E = 0.0;
    
    int proxyi ;
    Vector3d X;
    Vector3d N;
    Vector3i T;
    double e;
    
    for(int i=0; i<F.rows(); i++)
    {
        proxyi = R(i,0);
        if(l2orl21==0)
            e = l2error(clusters.at(proxyi),i);
        else
            e = l21error(clusters.at(proxyi),i);
        E+= e ;
//        cout<<"E: "<<e<<endl;
    }
    
    return E;
};

void bestfittingRegions(igl::opengl::glfw::Viewer &viewer){//},  double threshold){
//    int k= clusters.size();
    cout<<"going"<<"\n";
    vector<int> initialTriangles2;// = initialTriangles;
    initialTriangles2.clear();
    
    vector<Vector3d> clusterXi(k);
    vector<Vector3d> clusterNi(k);
    for(int i=0; i<k; i++){
        if(l2orl21==0)
        {
            clusterXi[i] = computeXi(clusters[i]);
            clusterNi[i] = computel2Ni(clusters[i]);
        }
        else
            clusterNi[i] = computel21Ni(clusters[i]);
    }
//    copy(initialTriangles.begin(), initialTriangles.end(), back_inserter(initialTriangles2));
    for (int i=0; i<initialTriangles.size(); i++)
            initialTriangles2.push_back(initialTriangles[i]);
    initialTriangles.clear();
    priority_queue<for_region, vector<for_region>, CompareError> Q;
    regTri = MatrixXi::Zero(F.rows(),2);
    VectorXi tris(k);
    VectorXd dist(k);
    double d;
    for(int i=0; i<k; i++)
        dist(i) = 9999.0;
    
        for(int i=0; i<F.rows(); i++)
        {
            int reg= regions(i,0);
//            cout<<" "<<reg<<" "<<clusters.at(reg).size()<<" \n";
//            clusterXi[reg] = computeXi(reg);
//            clusterNi[reg] = computeNi(reg);
            double errChk;
            if(l2orl21==0)
                errChk=l2error(clusterXi[reg], clusterNi[reg],  i);
            else
                errChk=l21error(clusterNi[reg],  i);
//            cout<<"cen: "<<computeXi(clusters.at(reg)).transpose()<<" \t nor: "<<computeNi(clusters.at(reg)).transpose()<<"\n";
            if(  errChk < dist(reg) )
            {
                dist(reg) = errChk;
                tris(reg) = i;
            }
        }
//        cout<<"d: "<<dist<<"\n"<<"tris"<<tris.transpose();
        clusters.clear();
        clusters.resize(k);
        regions = -MatrixXi::Ones(F.rows(),1);
        for(int i=0; i<k; i++ ){
            regions(tris(i))=i;
            initialTriangles.push_back(tris(i));
            clusters.at(i).push_back(tris(i));
            for(int r=0;r<neighbTri.at(tris(i)).size();r++){
                if(l2orl21==0)
                    Q.push(for_region( neighbTri.at(tris(i)).at(r), i, -l2error(clusterXi[i], clusterNi[i],  neighbTri.at(tris(i)).at(r) ) ));
                else
                    Q.push(for_region( neighbTri.at(tris(i)).at(r), i, -l21error(clusterNi[i],  neighbTri.at(tris(i)).at(r) ) ));
                
                    
            }
        }

        regTri = MatrixXi::Zero(F.rows(),2);
        while (Q.size()!=0) {
//            cout<<"going3"<<endl;
            for_region tr = Q.top();
            Q.pop();
            if(regTri(tr.face,0)==0 && regTri(tr.face,0)<3)
            {
                if(regions(tr.face)==-1)
                {
                    regions(tr.face) = tr.cluster;
                    regTri(tr.face,0)++;
                    regTri(tr.face,1)=tr.cluster;
                    clusters.at(tr.cluster).push_back(tr.face);
                    for(int j=0;j<neighbTri.at(tr.face).size();j++)
                    {
                        double asd;
                        if(l2orl21==0)
                            asd= l2error(clusterXi[tr.cluster], clusterNi[tr.cluster], neighbTri.at(tr.face).at(j));
                        else
                            asd= l21error(clusterNi[tr.cluster], neighbTri.at(tr.face).at(j));
                        //                    err+= asd;
//                        cout<<" \n"<<"face: "<<tr.face<<" "<<neighbTri.at(tr.face).at(j)<<" proxy: "<<tr.cluster<<" er: "<<asd<<endl;
                        Q.push(for_region( neighbTri.at(tr.face).at(j), tr.cluster, -asd ));
                    }
                }
            }
        }
        MatrixXd C;
        igl::jet(regions,true,C);
//        viewer.data().set_mesh(V, F);
//        viewer.append_mesh();
//        viewer.data().set_mesh(V, F);
        viewer.data().set_colors(C);
//        viewer.launch();
//    for(int i=0; i<initialTriangles.size(); i++)
//        cout<< initialTriangles.at(i)<<" \n";
    
    
    
//    return clusters;
};


void draw_ellipses(igl::opengl::glfw::Viewer& viewer, vector<vector<int>> clusters)
{
    
//    for(int i=0; i<F.rows(); i++)
//    {
//        for(int j=0; j<clusters.size(); j++)
//            if(vector_contains(clusters.at(j),i))
//                R(i,0) = j;
//    }
//
    int edges=20;
    MatrixXd ellipseV = MatrixXd::Zero(k*edges,3);
    MatrixXi ellipseF = MatrixXi::Zero(k*(edges-1),3);
    MatrixXi R = MatrixXi::Zero((edges-1)*k,3);
    for (int i=0; i<k; i++) {
        vector<Vector3d> r;
        for(int j=0; j<clusters.at(i).size(); j++)
        {
            for(int s=0; s<3; s++)
            {
//                cout<<clusters.at(i).size()<<"at ("<<i<<","<<j<<") = "<<clusters.at(i).at(j)<<endl;
                r.push_back((V.row(F(clusters.at(i).at(j),s))));
            }
        }
        MatrixXd M= MatrixXd::Zero(r.size(),3);
        for(int j=0; j<M.rows(); j++)
            M.row(j) = r.at(j);
        
        Vector3d mean = M.colwise().mean();
        for(int j=0; j<M.rows(); j++)
            M.row(j) -= mean;
        M= M.transpose()*M;
        EigenSolver<MatrixXd> eigen(M/r.size());
        MatrixXd Evectors = eigen.eigenvectors().real();
        MatrixXd Evalues = eigen.eigenvalues().real();
        
        Vector3d m1,m2;
        for(int j=0; j<3;j++)
        {
            if(Evalues(j) == Evalues.maxCoeff())
            {
                m1 = Evectors.col(j)*pow(Evalues(j),0.5);
                Evalues(j) = -Evalues(j);
                break;
            }
        }
        for(int j=0; j<3;j++)
        {
            if(Evalues(j) == Evalues.maxCoeff())
            {
                m2 = Evectors.col(j)*pow(Evalues(j),0.5);
            }
        }
        Vector3d normal;
        Vector3d Xi =computeXi(clusters.at(i)) ;
        
        if(l2orl21==0)
            normal= computel2Ni(clusters.at(i));
        else
            normal = computel21Ni(clusters.at(i));
        
        if(normal.dot(m1.cross(m2))<0)
        {
            m1 = -m1;
            ellipseV.row(edges*i) = Xi;
            for(int j=1; j<edges; j++)
            {
                double theta = j*2*PI / (edges-1);
                ellipseV.row(edges*i + j)= Xi+sin(theta)*m1 + cos(theta)*m2;
            }
        }
        else
        {
            ellipseV.row(edges*i) = Xi;
            for(int j=1; j<edges; j++)
            {
                double theta = j*2*PI / (edges-1);
                ellipseV.row(edges*i + j)= Xi+sin(theta)*m1 + cos(theta)*m2;
            }
        }
        
        for(int j=1; j<edges-1 ; j++)
        {
            ellipseF.row((edges-1)*i+j-1) = Vector3i(edges*i+j,edges*i,edges*i+j+1);
            R((edges-1)*i+j-1) = i;
        }
        ellipseF.row((edges-1)*(i+1)-1) = Vector3i(edges*i+edges-1,edges*i,edges*i+1);
        R((edges-1)*(i+1)-1) = i;
    }
    viewer.data().clear();
    viewer.data().set_mesh(ellipseV,ellipseF);
    MatrixXd C;
    
    igl::jet(R,true,C);
    viewer.data().set_colors(C);
}

int adding_anchors(HalfedgeDS he,double thresholdA,vector<vector<int>> clusters,int r,int anchor,VectorXi& anchors){
    int nA=0;
    MatrixXi R= MatrixXi::Zero(F.rows(),1);
    for(int i=0; i<F.rows(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(vector_contains(clusters.at(j),i))
            {
                R(i,0) = j;
            }
    }
    list<int> all_anchors;
    int edge = he.getOpposite(he.getEdge(anchor));
    while(R(he.getFace(he.getOpposite(edge)),0) == r || R(he.getFace(edge),0)!=r )
    {
        edge = he.getNext(he.getOpposite(edge));
//        cout<<"done"<<endl;
    }
    int fedge = edge;
    int nAnchor = he.getTarget(edge);
    while(nAnchor != anchor)
    {
        if(anchors(nAnchor)==1)
            all_anchors.push_back(nAnchor);
        edge = he.getNext(edge);
        while(R(he.getFace(he.getOpposite(edge)),0)==r)
        {
            edge = he.getNext(he.getOpposite(edge));
        }
        nAnchor = he.getTarget(edge);
    }
    all_anchors.push_back(anchor);
        
    edge = fedge;
    int x;
    while(!all_anchors.empty())
    {
        nAnchor = all_anchors.front();
        all_anchors.pop_front();
        int newV;// = he.getTarget(edge);
        double max=0.0;
        int maxV = -1;
        while(he.getTarget(edge)!= nAnchor)
        {
            edge = he.getNext(edge);
            while(R(he.getFace(he.getOpposite(edge)),0)==r)
            {
                edge = he.getNext(he.getOpposite(edge));
            }
            x= R(he.getFace(he.getOpposite(edge)),0);
            newV = he.getTarget(edge);
            double dist ;
            Vector3d nor1;
            Vector3d nor2;
            if(l2orl21 ==0)
            {
                nor1 = computel2Ni(clusters.at(r));
                nor2 = computel2Ni(clusters.at(x));
            }
            else{
                nor1 = computel21Ni(clusters.at(r));
                nor2 = computel21Ni(clusters.at(x));
            }
            
            double ang = (nor1.cross(nor2)).norm(); //norm or normalized?
            
            Vector3d s1 = V.row(newV) - V.row(anchor);
            Vector3d s2 = V.row(nAnchor) - V.row(anchor);
            
            if(s2.norm()==0)
                dist =0.0;
            else
            {
                double t = (s1.dot(s2))/s2.norm();
                dist = (s1 - t*s2/s2.norm()).norm()/s2.norm()*ang;
            }
            if(dist>max)
            {
                max = dist;
                maxV = newV;
            }
        }
        if(max > thresholdA)
        {
            anchors(maxV) = 1;
            nA++;
        }
        anchor = nAnchor;
    }
    
    return nA;
    
}

void draw_anchors(igl::opengl::glfw::Viewer& viewer, vector<vector<int>> clusters,HalfedgeDS he){
    
    
    MatrixXi R = MatrixXi::Zero(F.rows(),1);
    for(int i=0; i<F.rows(); i++)
    {
        for(int j=0; j<clusters.size(); j++)
            if(vector_contains(clusters.at(j),i))
                R(i,0) = j;
    }
    double thresholdA=0.4;
    VectorXi anchors;
    anchors.setZero(V.rows());
    vector<vector<int>> vertex_belongs_to(V.rows());
    for(int i=0; i<V.rows(); i++)
    {
//        for(int j=0; j<clusters.size(); j++)
//        {
//                if(vector_contains(clusters.at(j), i))
//                {
//                    vertex_belongs_to.at(i).push_back(j);
//                }
//        }
        int proxy,edge;
        edge = he.getEdge(i);
        proxy = R(he.getFace(edge),0);
        vertex_belongs_to.at(i).push_back(proxy);
        int pedge = he.getOpposite(he.getNext(edge));
        while (pedge!= edge) {
            proxy= R(he.getFace(pedge),0);
            if(!vector_contains(vertex_belongs_to.at(i),proxy))
               {
                vertex_belongs_to.at(i).push_back(proxy);
            }
               pedge = he.getOpposite(he.getNext(pedge));
        }
        if(vertex_belongs_to.at(i).size() > 2)
        {
            anchors(i) = 1;
        }
    }
    VectorXi visited;
    visited.setZero(k);
    int anc_size;
    for(int i=0; i<V.rows(); i++)
    {
        if(vertex_belongs_to.at(i).size() > 2)
        {
            for(int j=0; j<vertex_belongs_to.at(i).size(); j++ )
            {
                if(visited(vertex_belongs_to.at(i).at(j))==0)
                {
                    visited(vertex_belongs_to.at(i).at(j))=1;
                    anc_size= adding_anchors(he,thresholdA,clusters,vertex_belongs_to.at(i).at(j),i,anchors);
//                    while(anc_size>0)
//                    {
//                        anc_size = adding_anchors(he,thresholdA,clusters,vertex_belongs_to.at(i).at(j),i,anchors);
//                    }
                }
                
            }
        }
    }
    polygon_vertices.resize(k);
    for(int i=0; i<V.rows(); i++)
    {
        if(vertex_belongs_to.at(i).size()>2)
        {
            for(int j=0; j< vertex_belongs_to.at(i).size(); j++ )
            {
                int ss = vertex_belongs_to.at(i).at(j);
                if(!vector_contains(polygon_vertices.at(ss),i))
                {
                    vector<int> region_anchors;
                    region_anchors.push_back(i);
                    int fedge = he.getOpposite(he.getEdge(i));
                    while(R(he.getFace(he.getOpposite(fedge)),0) == ss || R(he.getFace(fedge),0)!=ss )
                    {
                        fedge = he.getNext(he.getOpposite(fedge));
                    }
                    int edge = fedge;
                    int new_anchor = he.getTarget(edge);
                    while(new_anchor!=i)
                    {
                        if(anchors(new_anchor)==1)
                        {
                            region_anchors.push_back(new_anchor);
                        }
                        edge = he.getNext(edge);
                        while(R(he.getFace(he.getOpposite(edge)),0)==ss)
                        {
                            edge = he.getNext(he.getOpposite(edge));
                        }
                        new_anchor = he.getTarget(edge);
                    }
                    polygon_vertices.at(ss).insert(polygon_vertices.at(ss).end(),region_anchors.begin(),region_anchors.end());
                }
            }
        }
    }
    for(int i=0; i<k; i++)
    {
        if(polygon_vertices.at(i).size()<3 && polygon_vertices.at(i).size()>0)
        {
            int edge = he.getNext(he.getEdge(polygon_vertices.at(i).at(0)));
            int new_anchor = he.getTarget(edge);
            polygon_vertices.at(i).push_back(new_anchor);
        }
    }
    for(int i=0; i<polygon_vertices.size(); i++)
    {
//        cout<<"at "<<i<<" the size is "<<polygon_vertices.at(i).size()<<endl;
        for(int j=0; j<polygon_vertices.at(i).size(); j++)
        {
//            cout<<"at "<<i<<" the anchor v's are "<<polygon_vertices.at(i)[j]<<endl;
            viewer.data().add_points(V.row(polygon_vertices.at(i).at(j)), Eigen::RowVector3d(1.,1,1));
        }
    }
    cout<<"done end"<<endl;
}

void triangulation(igl::opengl::glfw::Viewer& viewer)
{
    vector<Vector3d> trV;
    vector<Vector3i> trF;
    MatrixXd trVer;
    viewer.data().clear();
    for(int i=0; i<polygon_vertices.size(); i++)
    {
        Delaunay dt;
        vector<Point> pts;
        pts.clear();
        if(polygon_vertices.at(i).size()<3 && polygon_vertices.at(i).size()>0)
        {
            for(int j=0; j<polygon_vertices.at(i).size(); j++)
            {
                if(i<k-1)
                {
                    polygon_vertices.at(i+1).push_back(polygon_vertices.at(i).at(j));
                }
                else
                    polygon_vertices.at(i-1).push_back(polygon_vertices.at(i).at(j));
            }
//            polygon_vertices.at(i).clear();
        }
//        cout<<"at "<<i<<" the size is "<<polygon_vertices.at(i).size()<<endl;
        
            
            for(int j=0; j<polygon_vertices.at(i).size(); j++)
            {
                if(polygon_vertices.at(i).size()>2)
                {
                    dt.insert(Point(V(polygon_vertices.at(i).at(j),0), V(polygon_vertices.at(i).at(j),1), V(polygon_vertices.at(i).at(j),2)));
                    pts.push_back(Point(V(polygon_vertices.at(i).at(j),0), V(polygon_vertices.at(i).at(j),1), V(polygon_vertices.at(i).at(j),2)));
                }
                if( j< polygon_vertices.at(i).size()-1)
                {
                    MatrixXd nedge= MatrixXd::Zero(2,3);
                    nedge.row(0) = V.row(polygon_vertices.at(i).at(j));
                    nedge.row(1) = V.row(polygon_vertices.at(i).at(j+1));
//                    viewer.data(2).add_edges(nedge.row(0), nedge.row(1) , Eigen::RowVector3d(0,0,0));
                }
                else
                {
                    MatrixXd nedge= MatrixXd::Zero(2,3);
                    nedge.row(0) = V.row(polygon_vertices.at(i).at(j));
                    nedge.row(1) = V.row(polygon_vertices.at(i).at(0));
//                    viewer.data(2).add_edges(nedge.row(0), nedge.row(1) , Eigen::RowVector3d(0,0,0));
                }
                
                if( j== polygon_vertices.at(i).size()-1)
                    dt.insert_constraint(Point(V(polygon_vertices.at(i).at(j),0), V(polygon_vertices.at(i).at(j),1), V(polygon_vertices.at(i).at(j),2)), Point(V(polygon_vertices.at(i).at(0),0), V(polygon_vertices.at(i).at(0),1), V(polygon_vertices.at(i).at(0),2)));
                else
                    dt.insert_constraint(Point(V(polygon_vertices.at(i).at(j),0), V(polygon_vertices.at(i).at(j),1), V(polygon_vertices.at(i).at(j),2)), Point(V(polygon_vertices.at(i).at(j+1),0), V(polygon_vertices.at(i).at(j+1),1), V(polygon_vertices.at(i).at(j+1),2)));
                Polyhedron_3 convex_hull;
                CGAL::convex_hull_3(pts.begin(), pts.end(), convex_hull);
        //        std::cout << "The convex hull contains " << convex_hull.size_of_vertices() << " vertices" << std::endl;
                for (Polyhedron_3::Halfedge_iterator e = convex_hull.halfedges_begin(); e != convex_hull.halfedges_end(); e++)
                {
                    // retrieve source and target vertices of edge
                    auto source = e->vertex()->point();
                    auto target = e->opposite()->vertex()->point();
        //            cout<<" "<<Point(source).x()<<"\n";
                    MatrixXd nedge= MatrixXd::Zero(2,3);
                    nedge.row(0) = Vector3d(Point(source).x(), Point(source).y(), Point(source).z());
                    nedge.row(1) = Vector3d(Point(target).x(), Point(target).y(), Point(target).z());
        //            nedge.row(1) = Vector3d(Point(target));
        //            nedge.row(0) = Vector3d(source->point().x(), source->point().y(), source->point().z());
        //            nedge.row(1) = Vector3d(target->point().x(), target->point().y(), target->point().z());
                    
//                    viewer.data(2).add_edges(nedge.row(0), nedge.row(1) , Eigen::RowVector3d(0,0,0));
                }
                
                
            }
      
//        trV.clear();
        for (Delaunay::Vertex_iterator vit = dt.finite_vertices_begin();
             vit != dt.finite_vertices_end();
             ++vit)
        {
            Point p = vit->point();
            trV.push_back(Vector3d(p.x(), p.y(), p.z()));
        }
////        for(int j=0; j<trF.size();j++)
////        {
////            trFac.row(j)= trF.at(j);
////        }
        trVer=MatrixXd::Zero(trV.size(),3);
        for(int j=0; j<trV.size();j++)
        {
            trVer.row(j)= trV.at(j);
        }
        viewer.data(2).add_points(trVer, Eigen::RowVector3d(0.,0,0));
//
        Delaunay::Finite_edges_iterator eit;
        for (eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); eit++)
        {
            Delaunay::Vertex_handle source = eit->first->vertex((eit->second + 1) % 3);
            Delaunay::Vertex_handle target = eit->first->vertex((eit->second + 2) % 3);
////            auto source = eit->first->vertex(eit->second);
////                auto target = eit->first->vertex(eit->third);
////                    std::cout << "Edge: " << source->point() << " -> " << target->point() << std::endl;
            MatrixXd nedge= MatrixXd::Zero(2,3);
            nedge.row(0) = Vector3d(source->point().x(), source->point().y(), source->point().z());
            nedge.row(1) = Vector3d(target->point().x(), target->point().y(), target->point().z());
            viewer.data(2).add_edges(nedge.row(0), nedge.row(1) , Eigen::RowVector3d(0,0,0));
//           // cout << nedge.row(0)<<"\t"<<nedge.row(1)<<endl;
////                    viewer.data(1).add_edges(nedge.row(0), nedge.row(1) , Eigen::RowVector3d(0.5,0.5,0.5));
        }
    }
   
//    viewer.data().set_mesh(trVer, trFac);
        
        
        
    cout<<"done"<<endl;
    
}
double distance(Vector3d a , Vector3d b){
    return (b-a).norm();
    
}

MatrixXd dijkstra(MatrixXd G,int startnode) {
    int max = G.rows();
    int n= max;
   double cost[max][max],distance[max],pred[max];
   int visited[max],count,nextnode,i,j;
    double mindistance;
   for(i=0;i<n;i++)
       for(j=0;j<n;j++)
           if(G(i,j)==0)
               cost[i][j]=99999.0;
           else
               cost[i][j]=G(i,j);
    for(i=0;i<n;i++)
    {
        distance[i]=cost[startnode][i];
//      pred[i]=startnode;
        visited[i]=0;
   }
   distance[startnode]=0;
   visited[startnode]=1;
   count=1;
   while(count<n-1)
   {
       mindistance=99999.0;
       for(i=0;i<n;i++)
           if(distance[i]<mindistance&&!visited[i])
           {
               mindistance=distance[i];
               nextnode=i;
           }
       visited[nextnode]=1;
       for(i=0;i<n;i++)
           if(!visited[i])
               if(mindistance+cost[nextnode][i]<distance[i])
               {
                   distance[i]=mindistance+cost[nextnode][i];
//         pred[i]=nextnode;
               }
       count++;
   }
    MatrixXd targets = MatrixXd::Zero(n,1);
    for(i=0;i<n;i++)
    {
        targets(i,0) = distance[i];
    }
//      cout<<"
//Path="<<i;
//      j=i;
//      do {
//         j=pred[j];
//         cout<<"<-"<<j;
//      }while(j!=startnode);
//   }
    return targets;
}



void triangulation2(igl::opengl::glfw::Viewer& viewer)
{
    viewer.data().clear();
//    MatrixXd newVer= MatrixXd::Zero(V.rows()*50,3);
//    MatrixXi newFac= MatrixXi::Zero(F.rows()*10,3);
    int ind = 0, facind=0;
    for(int i=0; i<clusters.size(); i++)
//    for(int i=0; i<3; i++)
    {
        MatrixXd Vofechcluster;
        
        vector<int> verts;
        set<int> target;
        verts.clear();
        for(int j=0; j<clusters[i].size(); j++)
        {
            if(!vector_contains(verts,  F(clusters.at(i).at(j),0))   )
                verts.push_back(F(clusters.at(i).at(j),0));
            if(!vector_contains(verts,F(clusters.at(i).at(j),1)))
                verts.push_back(F(clusters.at(i).at(j),1));
            if(!vector_contains(verts,F(clusters.at(i).at(j),2)))
                verts.push_back(F(clusters.at(i).at(j),2));
            target.insert(F(clusters.at(i).at(j),0));
            target.insert(F(clusters.at(i).at(j),1));
            target.insert(F(clusters.at(i).at(j),2));
        }
        
        
        MatrixXd anch_weights = MatrixXd::Zero(polygon_vertices[i].size(),polygon_vertices[i].size());
        for(int j=0; j<anch_weights.rows(); j++)
        {
//            for(int jj=0; jj<anch_weights.cols(); jj++)
            if(j==anch_weights.rows()-1)
            {
                anch_weights(j,0) = distance(V.row(polygon_vertices[i][j]) , V.row(polygon_vertices[i][0]));
                anch_weights(0,j) = distance(V.row(polygon_vertices[i][0]) , V.row(polygon_vertices[i][j]));
            }
            else
            {
                anch_weights(j,j+1) = distance(V.row(polygon_vertices[i][j]) , V.row(polygon_vertices[i][j+1]));
                anch_weights(j+1,j) = distance(V.row(polygon_vertices[i][j+1]) , V.row(polygon_vertices[i][j]));
            }
        }
        MatrixXd anchor_weights = MatrixXd::Zero(polygon_vertices[i].size(),polygon_vertices[i].size());
        for(int j=0; j<polygon_vertices[i].size(); j++)
        {
            anchor_weights.col(j)=dijkstra(anch_weights, j);
        }
//        cout<<anch_weights<<"\nxxxxxxxxxxxxxxx\n"<<endl;
//        cout<<anchor_weights<<endl;

        
        vector<int> temps(polygon_vertices[i].size());
        for(int j=0; j<polygon_vertices[i].size() ; j++)
            temps.push_back(j);
        
        
//        while(temps.size()>3)
//        {
//            double minimum = 9999999.0;
//            int minimum_index;
//            for( int j =0 ; j< temps.size(); j++)
//            {
//                if(j==0)
//                {
//                    if(anchor_weights(temps[temps.size()-1],temps[j+1]) <minimum )
//                    {
//                        minimum = anchor_weights(temps[temps.size()-1],temps[j+1]);
//                        minimum_index = j;
//                    }
//                }
//                else if(j==temps.size()-1)
//                {
//                    if(anchor_weights(temps[j-1],temps[0]) <minimum )
//                    {
//                        minimum = anchor_weights(temps[j-1],temps[0]);
//                        minimum_index = j;
//                    }
//                }
//                else{
//                    if(anchor_weights(temps[j-1],temps[j+1]) <minimum )
//                    {
//                        minimum = anchor_weights(temps[j-1],temps[j+1]);
//                        minimum_index = j;
//                    }
//
//                }
////                if()
////                cout<<"done"<<endl;
//            }
//            if(minimum_index==0)
//            {
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[temps.size()-1]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index+1]]);
//                newFac.row(facind) = Vector3i(ind-2, ind-1,ind);
//                ind++;
//                facind++;
//            }
//            else if(minimum_index == temps.size()-1)
//            {
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index-1]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[0]]);
//                newFac.row(facind) = Vector3i(ind-2, ind-1,ind);
//                ind++;
//                facind++;
//            }
//            else
//            {
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index-1]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index]]);
//                ind++;
//                newVer.row(ind) = V.row(polygon_vertices[i][temps[minimum_index+1]]);
//                newFac.row(facind) = Vector3i(ind-2, ind-1,ind);
//                ind++;
//                facind++;
//            }
//            vector<int>::iterator it;
//            it = temps.begin()+minimum_index;
//            temps.erase(it);
//
//        }
//        newVer.row(ind) = V.row(polygon_vertices[i][temps[0]]);
//        ind++;
//        newVer.row(ind) = V.row(polygon_vertices[i][temps[1]]);
//        ind++;
//        newVer.row(ind) = V.row(polygon_vertices[i][temps[2]]);
//        newFac.row(facind) = Vector3i(ind-2, ind-1,ind);
//        ind++;
//        facind++;
//
        
        
        
        for(int j=0; j< polygon_vertices[i].size(); j++)
        {
//            int edge = he.getEdge(polygon_vertices[i][j]);
//            while(vector_contains(polygon_vertices[i],he.getTarget(he.getOpposite(edge)) ))
//            {
//                edge = he.getOpposite(he.getNext(edge));
//            }
//            if(!vector_contains(verts, he.getTarget(he.getOpposite(edge)  ) ) )
//            {
//                bedges.push_back(edge);
//                viewer.data().add_edges(V.row(he.getTarget(edge)))
//            }
            if(j== polygon_vertices[i].size()-1 )
            viewer.data().add_edges(V.row(polygon_vertices[i][j]), V.row(polygon_vertices[i][0]), Eigen::RowVector3d(0,1.0,1.0));
            else
            viewer.data().add_edges(V.row(polygon_vertices[i][j]), V.row(polygon_vertices[i][j+1]), Eigen::RowVector3d(0,1.0,1.0));
                
        }
        
        
        
        
//        verts.assign(target.begin(),target.end());
//        cout<<"size: "<<target.size()<<" "<<verts.size()<<endl;
//        for(int xx : target)
//            cout<<"targ "<<xx<<" \n";
//        for(int xx : verts)
//            cout<<"vert "<<xx<<" \n";
        Vofechcluster =MatrixXd::Zero(verts.size(),3);
        
        MatrixXi connections = MatrixXi::Zero(verts.size(),verts.size());
        MatrixXd connections_weight = MatrixXd::Zero(verts.size(),verts.size());
        
        set<int> targets;//(0,Vofechcluster.rows());//(verts.begin(),verts.end());
        for(int j=0; j< verts.size(); j++)
        {
            targets.insert(j);
            Vofechcluster.row(j) = V.row(verts.at(j));
        }
        
        MatrixXi anchorind;
        anchorind = MatrixXi::Zero(polygon_vertices[i].size(),1);
        for(int j=0; j<polygon_vertices[i].size(); j++)
        {
            for(int it=0; it<Vofechcluster.rows(); it++)
            {
                if(V.row(polygon_vertices[i].at(j)) ==  Vofechcluster.row(it))
                {
                    anchorind(j,0) = it;
//                    cout<< polygon_vertices[i].at(j) <<" " <<Vofechcluster(it,0)<< " " <<anchorind(j,0)<<endl;
                }
            }
        }
//        cout<<"xxxxxxxxxxxx\n"<<Vofechcluster;
        MatrixXd flooding;
        flooding=-MatrixXd::Ones(verts.size(),2);
        flooding.col(1) = -999.0 * flooding.col(1);
        for(int j=0; j<clusters[i].size(); j++)
        {
            int x,y,z;
            //            x= V.row(F(clusters[i][j],0));
            //            y= V.row(F(clusters[i][j],1));
            //            z= V.row(F(clusters[i][j],2));
            for(int it=0; it<Vofechcluster.rows(); it++)
            {
                if(V.row(F(clusters[i][j],0)) == Vofechcluster.row(it))
                {
                    x=it;
                }
                if(V.row(F(clusters[i][j],1)) == Vofechcluster.row(it))
                {
                    y=it;
                }
                if(V.row(F(clusters[i][j],2)) == Vofechcluster.row(it))
                {
                    z=it;
                }
            }
            connections(x,y) = 1;
            connections(y,x) = 1;
            connections(x,z) = 1;
            connections(z,x) = 1;
            connections(y,z) = 1;
            connections(z,y) = 1;
            
            connections_weight(x,y) = distance(Vofechcluster.row(x), Vofechcluster.row(y));
            connections_weight(y,x) = distance(Vofechcluster.row(x), Vofechcluster.row(y));
            connections_weight(x,z) = distance(Vofechcluster.row(x), Vofechcluster.row(z));
            connections_weight(z,x) = distance(Vofechcluster.row(x), Vofechcluster.row(z));
            connections_weight(y,z) = distance(Vofechcluster.row(z), Vofechcluster.row(y));
            connections_weight(z,y) = distance(Vofechcluster.row(z), Vofechcluster.row(y));
            
            
        }
//        cout << "\n\n"<<"\n"<<connections;

        vector<vector<int>> VV(verts.size());
        vector<double> weights;
        for(int ij=0; ij<connections.rows(); ij++)
        {
            for(int jj=0; jj<connections.cols(); jj++ )
            {
                if(connections(ij,jj) == 1)
                {
                    VV[ij].push_back(jj);
//                    VV[jj].push_back(ij);
                    weights.push_back(connections_weight(ij,jj));
                }
            }
        }
        
        for(int ij=0; ij<anchorind.rows(); ij++)
        {
//            set<int> target = targets;
//            target.erase(anchorind(ij,0));
            Eigen::MatrixXd min_distance;
//                min_distance.setZero(targets.size());
              Eigen::VectorXi previous;
            min_distance = dijkstra(connections_weight, anchorind(ij,0));
//              igl::dijkstra(anchorind(ij,0), targets, VV, weights, min_distance, previous);
//            cout<<anchorind(ij,0)<<"gg"<<min_distance.size()<<"\t"<<target.size()<<"\n"<<min_distance<<endl;
//            cout<<"at ij= "<<ij<<" anchorpoint is: "<<anchorind(ij,0)<<" no of points in the region are: "<<verts.size() <<endl;
//                double md=0.0;
//            cout<<"at ij= "<<ij<<" anchorpoint is: "<<anchorind(ij,0)<<"\n......\n";
//            for(int xx : targets)cout<<xx<<" "<<endl;
            for(int xx=0; xx<targets.size(); xx++)
            {
//                cout<<"going"<<endl;
                if(min_distance(xx)< flooding(xx,1))
                {
//                    cout<<"\ndist: "<<min_distance(xx)<<" floodx1: "<<flooding(xx,1)<<" from "<<anchorind(ij,0)<<" to: "<<verts.at(xx)<<endl;
//                    for(int ss=0; ss<VV[anchorind(ij,0)].size(); ss++)cout<<VV[anchorind(ij,0)][ss]<<" ";
                    flooding(xx,0)=anchorind(ij,0);
                    flooding(xx,1)=min_distance(xx);
                }

            }
//            cout << "\n\n "<<anchorind(ij,0)<<" "<<verts.size()<<"\n"<< flooding;
        }
//        vector<vector<int>> buildfaces;
//        buildfaces.clear();
//        buildfaces.resize(clusters[i].size());
        
        for(int j=0; j< clusters[i].size(); j++ )
        {
            int xx,yy,zz;
            for(int jj=0; jj<Vofechcluster.rows(); jj++)
            {
                
                if(V.row( F(clusters[i][j],0)  ) == Vofechcluster.row(jj))
                    xx=jj;
                if(V.row( F(clusters[i][j],1)  ) == Vofechcluster.row(jj))
                    yy=jj;
                if(V.row( F(clusters[i][j],2)  ) == Vofechcluster.row(jj))
                    zz=jj;
                //                buildfaces[j].push_back(xx);
                //                buildfaces[j].push_back(yy);
                //                buildfaces[j].push_back(zz);
            }
            if(flooding(xx,0)!=flooding(yy,0) && flooding(xx,0)!= flooding(zz,0)&& flooding(yy,0)!= flooding(zz,0))
            {
                //                                cout<<"going"<<endl;
//                viewer.data(2).add_points(Vofechcluster.row(xx), Eigen::RowVector3d(0.,0,1.0) );
//                viewer.data(2).add_points(Vofechcluster.row(yy), Eigen::RowVector3d(0.,0,1.0) );
//                viewer.data(2).add_points(Vofechcluster.row(zz), Eigen::RowVector3d(0.,0,1.0) );
//                viewer.data().add_edges(Vofechcluster.row(xx),Vofechcluster.row((yy)), Eigen::RowVector3d(1.,0,1.0) );
//                viewer.data().add_edges(Vofechcluster.row(yy),Vofechcluster.row((zz)), Eigen::RowVector3d(1.,0,1.0) );
//                viewer.data().add_edges(Vofechcluster.row(zz),Vofechcluster.row((xx)), Eigen::RowVector3d(1.,0,1.0) );
//                viewer.data().add_edges(Vofechcluster.row(xx),Vofechcluster.row(flooding(xx,0)), Eigen::RowVector3d(0.,1,1) );
//                viewer.data().add_edges(Vofechcluster.row(yy),Vofechcluster.row(flooding(yy,0)), Eigen::RowVector3d(0.,1,1) );
//                viewer.data().add_edges(Vofechcluster.row(zz),Vofechcluster.row(flooding(zz,0)), Eigen::RowVector3d(0.,1,1) );
                viewer.data().add_edges(Vofechcluster.row(flooding(yy,0)),Vofechcluster.row(flooding(xx,0)), Eigen::RowVector3d(0.,1,1) );
                viewer.data().add_edges(Vofechcluster.row(flooding(zz,0)),Vofechcluster.row(flooding(yy,0)), Eigen::RowVector3d(0.,1,1) );
                viewer.data().add_edges(Vofechcluster.row(flooding(xx,0)),Vofechcluster.row(flooding(zz,0)), Eigen::RowVector3d(0.,1,1) );
//                newVer.row(ind) = Vofechcluster.row(flooding(xx,0));
//                ind++;
//                newVer.row(ind) = Vofechcluster.row(flooding(yy,0));
//                ind++;
//                newVer.row(ind) = Vofechcluster.row(flooding(zz,0));
//
//                newFac.row(facind) = Vector3i(ind-2, ind-1,ind);
//                ind++;
//                facind++;
            }
            //                    (flooding(jj,0))))
            //                F(clusters[i][j],1)
            //                F(clusters[i][j],2)
            
        }


        
        
        
        
        
        
    }
    cout<<"donedone"<<endl;
//    newFac = newVer.rowwise().filter([](const Eigen::RowVector3d& row) {
//      return (row.array() != 0).count() == 3;
//    });
//    ACIAIC
//    cout<<"size of vertices: "<<newVer.rows()<<endl;
//    MatrixXd newVer2;
//    MatrixXi IA,IC;
//    igl::unique_rows(newVer, newVer2, IA,IC);
//    cout<<" vertices: "<<newVer<<endl;

//    viewer.data().clear();
//    cout<<"done"<<endl;
//    viewer.data().set_mesh(newVer, newFac);
    cout<<"done2"<<endl;
//    viewer.data(2).set_mesh(newVer, newFac);
//    Eigen::MatrixXd V(4,3);
//    V << 0,0,0, 1,0,0, 2,0,0, 3,0,0;
//    std::vector<std::vector<int> > VV(4);
//    VV[0] = {1};
//    VV[1] = {2};
//    VV[2] = {3};
//    VV[3] = {0};
//    Eigen::VectorXd weights(4);
//    weights << 1,2,1,2;
//
//    // Define the starting point of the shortest path search
//    int source = 0;
//
//    // Define the target vertices for the shortest path search
//    std::set<int> targets = {1,2,3};
//    MatrixXd anc(polygon_vertices[2].size(),3);
//    set<int> anchset;
////    vector<int> forvv(polygon_vertices[2].size());
////    forvv.setZero();
////    MatrixXi forvv = MatrixXi::Zero(clusters[2].size(),clusters[2].size());
//    vector<vector<int>> VV(clusters[2].size());
//    MatrixXd testing2(clusters[2].size()*3,3);
//    VectorXd mindist;
//    VectorXi prev;
//
//    for(int i=0; i< polygon_vertices[2].size(); i++)
//    {
//        anc.row(i) = V.row(polygon_vertices[2].at(i));
//        anchset.insert(polygon_vertices[2].at(i));
//    }
//    int j=0;
//    for(int i=0; i<clusters[2].size(); i++)
//    {
//        testing2.row(j) = V.row(F(clusters[2][i],0));
//        testing2.row(j+1) = V.row(F(clusters[2][i],1));
//        testing2.row(j+2) = V.row(F(clusters[2][i],2));
//        j++;
//        VV(F(clusters[2][i],0))
//    }
    
//    igl::dijkstra(testing2,VV,polygon_vertices[2].at(0),anchset, mindist , prev );
//    igl::dijkstra(V,);
//    cout<<mindist;
//    MatrixXi flooding =MatrixXi::Zero(V.rows(),2);
//    for(int i=0; i<clusters.size(); i++)
//    {
//
//        vector<int> verts;
//        verts.clear();
//        for(int j=0; j<clusters.at(i).size(); j++)
//        {
//            if(!vector_contains(verts,F(clusters.at(i).at(j),0)))
//            verts.push_back(F(clusters.at(i).at(j),0));
//            if(!vector_contains(verts,F(clusters.at(i).at(j),1)))
//            verts.push_back(F(clusters.at(i).at(j),1));
//            if(!vector_contains(verts,F(clusters.at(i).at(j),2)))
//            verts.push_back(F(clusters.at(i).at(j),2));
//
//        }
//        for(int ii=0; ii<vert.size(); ii++)
//        {
//            igl::dijkstra();
//            flooding(ii , 0) =
//            flooding(F(clusters.at(i).at(j),0) , 1) =
//
//        }
//    }
}


int iterations =0;
double error = 100. , prior_error = 999.0;

//errs.push_back(error);
// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    switch(key){
        case '1':
            bestfittingRegions( viewer);//, threshold);
            iterations++;
            error = disrtortionError(clusters);
            errs.push_back(error);
            cout<<"Error at "<<iterations<< "th iteration is: "<<error<<"\n";
//            viewer.data().set_normals(N_faces);
            return true;
        case '2':
        {
            
            while (iterations <100 && fabs(error - prior_error)> threshold ) {
                
                prior_error = error;
                iterations++;
               
                if (vector_contains(errs,error) && errs.size()>6) {
                    cout<<"error is repeating\n";
                    break;
                }
                bestfittingRegions( viewer);
                errs.push_back(error);
                error = disrtortionError(clusters);
                
                cout<<"Error at "<<iterations<< "th iteration is: "<<error<<"\n";
            }
            
            
        }
            return true;
        case '3':
            draw_ellipses(viewer,clusters);
            return true;
        case '4':
            draw_anchors(viewer,clusters,*he);
            return true;
        case '5':
            triangulation2(viewer);
            return true;
        case '6':
            viewer.data().set_mesh(V,F);
            return true;
//        case '7':
//            viewer.data().set_normals(lib_N_vertices);
//            return true;
//        case '8':
//            viewer.data().set_normals(he_N_vertices);
//            return true;
        default: break;
    }
    return false;
}


// ------------ main program ----------------
int main(int argc, char *argv[]) {

//	if(argc<2) {
//		std::cout << "Error: input file required (.OFF)" << std::endl;
//		return 0;
//	}
//	std::cout << "reading input file: " << argv[1] << std::endl;

//	igl::readOFF(argv[1], V, F);
    k=160;
    threshold = 0.001;
    l2orl21=1;
    string file = "../data/bunny.off";
    if (argc>=2) {
      string w = argv[1];
      file = "../data/" + w + ".off";
    }
    if (argc>=3) {
      k = atoi(argv[2]);
    }
    if (argc>=4) {
      l2orl21 = atoi(argv[3]);
    }
    if (argc>=5) {
      threshold = atof(argv[4]);
    }
    
//        igl::readOFF("../data/homer1.off",V,F);
    igl::readOFF(file, V, F);
    
	//print the number of mesh elements
        //std::cout << "Points: " << V.rows() << std::endl;


       HalfedgeBuilder* builder=new HalfedgeBuilder();  //

       HalfedgeDS he2=builder->createMesh(V.rows(), F);  //
       
    he = &he2;
    
    //initial seeding
//    vector<int,double> forcomparison; // face, its cluster/region , l2 error
    priority_queue<for_region, vector<for_region>, CompareError> Q;
    
    SparseMatrix<int> a(F.rows(),F.rows());
//    vector<vector<int>> neighbTri;
    igl::facet_adjacency_matrix(F,a);
//    neighbTri.resize(a.outerSize());
    for (int i=0; i<a.outerSize(); ++i){
        vector<int> f;
        for (SparseMatrix<int>::InnerIterator it(a,i); it; ++it)
        {
            //        it.value();
            f.push_back(it.row());   // row index
            
            //        it.col();   // col index (here it is equal to i)
            //        it.index(); // inner index, here it is equal to it.row()
//            cout<<i<<" "<<it.value()<<" "<<it.row()<<" "<<it.col()<<" "<<it.index()<<"\n";
        }
        neighbTri.push_back(f);
        
    }
    
    
    regTri = MatrixXi::Zero(F.rows(),2); //1st col for marking&no of entries 2nd col for cluster id //for checking proxy assignment
    Vector3d clusterXi[k];
    Vector3d clusterNi[k];
//        for(int i=0; i<neighbTri.size(); i++)
//        {
////            for(int j=0; j<neighbTri.at(i).size();j++)
////            {
//            cout<<"at i= "<<i<<"size: "<<neighbTri.at(i).size()<<endl;;
////            }
//        }
    for(int i=0;i<k;i++)
    {
        int ran= rand() % F.rows();
        
//        cout<<"going "<<i<<"\n";
        if(!vector_contains(initialTriangles, ran))
        {
            initialTriangles.push_back(ran);
            vector<int> triangles;
            triangles.push_back(ran);
//            clusters.push_back(triangles);
//            cout<<"going "<<ran<<"\n";
            if(l2orl21==0)
            {
                clusterXi[i] = computeXi(ran);
                clusterNi[i] = computel2Ni(ran);
            }
            else
                clusterNi[i] = computel2Ni(ran);
            
            for(int j=0; j<neighbTri.at(ran).size(); j++)
            {
//                cout<<"cent "<<clusterXi[i].transpose()<<"\t"<<"nor "<< clusterNi[i].transpose()<<"\n";
//                cout<<"before "<<i<<" "<<ran<<" "<<neighbTri.at(ran).at(j)<<endl;
                double l;
                if(l2orl21==0)
                    l=l2error( clusterXi[i],clusterNi[i],neighbTri.at(ran).at(j));
                else
                    l=l21error(clusterNi[i],neighbTri.at(ran).at(j));
//                cout<<"error: "<<l<<"\n";
                Q.push(for_region( neighbTri.at(ran).at(j), i, -l  ));
//                cout<<"after "<<i<<" "<<ran<<" "<<neighbTri.at(ran).at(j)<<endl;
            }
//            initialTriangles2.push_back(neighbTri.at(ran).at(0));
//            Q.push(for_region( neighbTri.at(ran).at(0), i, l2error(neighbTri.at(ran).at(0)) ));
//            initialTriangles2.push_back(neighbTri.at(ran).at(1));
//            Q.push(for_region( neighbTri.at(ran).at(1), i, l2error(neighbTri.at(ran).at(1)) ));
//            initialTriangles2.push_back(neighbTri.at(ran).at(2));
//            Q.push(for_region( neighbTri.at(ran).at(2), i, l2error(neighbTri.at(ran).at(2)) ));
            
        }
        else
        {
            i--;
        }
    }
    clusters.resize(k);
//    cout<<"going"<<endl;
//    for(int i=0; i<initialTriangles.size(); i++)
//        cout<< initialTriangles.at(i) <<" \n";
    
    while (Q.size()!=0) {
       for_region tr = Q.top();
//        cout<<"face: "<<tr.face<<"\n";
        Q.pop();
        if(regTri(tr.face,0)==0 && regTri(tr.face,0)<3)
        {
            regTri(tr.face,0)++;
            regTri(tr.face,1)=tr.cluster;
//            cout<<"going "<<tr.cluster<<" "<<clusters.size()<<"\n";
            clusters.at(tr.cluster).push_back(tr.face);
            for(int j=0;j<neighbTri.at(tr.face).size();j++)
            {
                double l;
                if(l2orl21==0)
                    l=l2error( clusterXi[tr.cluster],clusterNi[tr.cluster],neighbTri.at(tr.face).at(j));
                else
                    l=l21error(clusterNi[tr.cluster],neighbTri.at(tr.face).at(j));
                
                Q.push(for_region( neighbTri.at(tr.face).at(j), tr.cluster , -l ));
                
            }
//            Q.push(for_region( neighbTri.at(tr.face).at(0), tr.cluster, l2error(clusters.at(tr.cluster), neighbTri.at(tr.face).at(0)) ));
//            Q.push(for_region( neighbTri.at(tr.face).at(1), tr.cluster, l2error(clusters.at(tr.cluster), neighbTri.at(tr.face).at(1)) ));
//            Q.push(for_region( neighbTri.at(tr.face).at(2), tr.cluster, l2error(clusters.at(tr.cluster), neighbTri.at(tr.face).at(2)) ));
        }
    }
    

//    cout<<he.sizeOfVertices()<<" "<< he.sizeOfFaces() <<" \n";
//    for (int x=0; x<clusters.size(); x++) {
////        for (int y=0; y<clusters.at(3).size(); y++) {
////            cout<<"3"<<" "<<clusters.at(3).size()<<" "<<clusters.at(3).at(y)<<endl;
//        cout<<x<<" "<<clusters.at(x).size()<<endl;
////        }
//    }
    cout<<"Error: "<<disrtortionError(clusters)<<"\n";
    
//    he.print();
    
    

	// Compute per-face normals
//	igl::per_face_normals(V,F,N_faces);


	// Compute per-vertex normals
//	igl::per_vertex_normals(V,F,N_vertices);

 ///////////////////////////////////////////

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer

  viewer.callback_key_down = &key_down;
//  viewer.data().show_lines = false;
  viewer.data().set_mesh(V, F);  //
//  viewer.data().set_normals(N_faces);  //
  std::cout<<
    "Press '1' To perform one iteration."<<std::endl<<
	"Press '2' To perform until convergence."<<std::endl<<
    "Press '3' To generate proxies ."<<std::endl<<
    "Press '4' To create anchor vertices"<<std::endl<<
	"Press '5' To create Final output"<<std::endl;

    regions = -MatrixXi::Ones(F.rows(),1);
    for(int i=0; i<F.rows(); i++)
    {
        for(int j=0; j<clusters.size();j++)
        {
            if(vector_contains(clusters.at(j),i)){
                regions(i,0)=j;
                //        cout<<neighbTri.at(i).at(j)<<" ";
                //        cout<<"\n";
            }
            
        }
    }
    VectorXd Z;
    Z.setZero(V.rows(),1);

    Z.col(0)=V.col(1);

    MatrixXd C;
  // Assign per-vertex colors
	igl::jet(Z,true,C);
	viewer.data().set_colors(C);  // Add per-vertex colors
  //viewer.core(0).align_camera_center(V, F);  //not needed
//    viewer.append_mesh();
	viewer.launch(); // run the viewer
    
//    clusters = bestfittingRegions( viewer, clusters, he , threshold, regTri);
}
