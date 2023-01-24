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

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "HalfedgeBuilder.cpp"



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
int k=180;
double threshold = 0.0000000001;//0.000000001;//0.001;

vector<vector<int>> clusters; // no of proxied regions or clusters
vector<int> initialTriangles;
vector<double> errs;



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
    return (  (v2-v1).cross(v3-v1) / (  (v2-v1).cross(v3-v1).norm() )  );
}

Vector3d triangleNormal(MatrixXd t){
    Vector3d v1 = t.row(0);
    Vector3d v2 = t.row(1);
    Vector3d v3 = t.row(2);
    return (  (v2-v1).cross(v3-v1) / (  (v2-v1).cross(v3-v1).norm() )  );
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
Vector3d computeNi(vector<MatrixXd> t){
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
Vector3d computeNi(MatrixXd t){
    
    return triangleNormal(t);
}

Vector3d computeNi(int t){
    MatrixXd tr=MatrixXd::Zero(3,3);
    for( int i=0;i<3;i++)
    tr.row(i)= V.row(F(t,i));
    return triangleNormal(tr);
}

Vector3d computeNi(vector<int> t)
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
    return computeNi(tr);
}
double orthoDistance(Vector3d Xi, Vector3d Ni, Vector3d v){
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
    return l2error(computeXi(tr) , computeNi(tr) , t);
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
    
    return l2error(computeXi(t), computeNi(t), t);
}



double l21error(Vector3d Ni, MatrixXd T){
    return ((triangleNormal(T) - Ni).squaredNorm()) * areaOfTriangle(T);
}

double l21error(int f){
    
    MatrixXd t = MatrixXd::Zero(3,3);
    t.row(0)= V.row(F(f,0));
    t.row(1)= V.row(F(f,1));
    t.row(2)= V.row(F(f,2));
    
    return l21error( computeNi(t), t);
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
    return l21error(computeNi(tr) , t);
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

    int p = clusters.size();
    int f = F.rows();
    MatrixXi R(f,1);
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
    
    for(int i=0; i<f; i++)
    {
        proxyi = R(i,0);
        e = l2error(clusters.at(proxyi),i);
        E+= e ;
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
        clusterXi[i] = computeXi(clusters[i]);
        clusterNi[i] = computeNi(clusters[i]);
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
    
//    while(error > threshold)
     // threshold = 0.000000000001
//        cout<<"going"<<endl;
        //        err=0;
//        cout<<"regions xxxxxxxx\n"<< regions<<"\nregions xxxxxxxx\n";
        
        for(int i=0; i<F.rows(); i++)
        {
            int reg= regions(i,0);
//            cout<<" "<<reg<<" "<<clusters.at(reg).size()<<" \n";
//            clusterXi[reg] = computeXi(reg);
//            clusterNi[reg] = computeNi(reg);
            double errChk=l2error(clusterXi[reg], clusterNi[reg],  i);
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
                Q.push(for_region( neighbTri.at(tris(i)).at(r), i, -l2error(clusterXi[i], clusterNi[i],  neighbTri.at(tris(i)).at(r) ) ));
            }
        }

//            for(int i=0; i<initialTriangles.size(); i++)
//                cout<<"init tri: "<< initialTriangles.at(i) <<" \n";
//        cout<<"going2"<<endl;
        regTri = MatrixXi::Zero(F.rows(),2);
//        cout<<"error= "<<err<<"\n";
//        if(error==err)
//            cont++;
//        if(cont==1)
//            err=0;
//        error=err;
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
                        double asd= l2error(clusterXi[tr.cluster], clusterNi[tr.cluster], neighbTri.at(tr.face).at(j));
                        //                    err+= asd;
//                        cout<<" \n"<<"face: "<<tr.face<<" "<<neighbTri.at(tr.face).at(j)<<" proxy: "<<tr.cluster<<" er: "<<asd<<endl;
                        Q.push(for_region( neighbTri.at(tr.face).at(j), tr.cluster, -asd ));
                    }
                }
            }
            
//            if(regTri(tr.face,0)<3 && regions(tr.face)==-1)
//            {
////             if(regions(tr.face)==-1)
////             {
//                 regions(tr.face) = tr.cluster;
//                regTri(tr.face,0)++;
//                regTri(tr.face,1)=tr.cluster;
//                clusters.at(tr.cluster).push_back(tr.face);
//                for(int j=0;j<3;j++)
//                {
//                    double asd= l2error(computeXi(clusters.at(tr.cluster)), computeNi(clusters.at(tr.cluster)), neighbTri.at(tr.face).at(j));
////                    err+= asd;
//                    Q.push(for_region( neighbTri.at(tr.face).at(j), tr.cluster, asd ));
//                }
//            }
           
        
        }
        
        
//        for (int x=0; x<clusters.size(); x++) {
////            for (int y=0; y<clusters.at(x).size(); y++) {
////                cout<<x<<" "<<clusters.at(x).size()<<" "<<clusters.at(x).at(y)<<endl;
//            cout<<x<<" "<<clusters.at(x).size()<<endl;
////            }
//        }
//        cout<<"Error: "<<disrtortionError(clusters)<<"\n";
//        VectorXd Z;
//
//        Z.setZero(F.rows(),1);
//        for(int i=0; i<F.rows(); i++)
//        {
//            for(int j=0; j<clusters.size(); j++)
//            if(vector_contains(clusters.at(j),i))
//            {
//                Z(i,0)=j;
//            }
//        }
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
            
            while (iterations <100 && fabs(error - prior_error)> 0.0001 ) {
                if (vector_contains(errs,error)==true) {
                    cout<<"error is repeating\n";
                    break;
                }
                prior_error = error;
                iterations++;
                bestfittingRegions( viewer);
                errs.push_back(error);
                error = disrtortionError(clusters);
                
                cout<<"Error at "<<iterations<< "th iteration is: "<<error<<"\n";
            }
            
            
        }
            return true;
        case '3':
            viewer.data().set_normals(N_vertices);
            return true;
        case '4':
            viewer.data().set_normals(lib_N_vertices);
            return true;
        case '5':
            viewer.data().set_normals(he_N_vertices);
            return true;
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
    
        igl::readOFF("../data/nefertiti.off",V,F);

	//print the number of mesh elements
        //std::cout << "Points: " << V.rows() << std::endl;


       HalfedgeBuilder* builder=new HalfedgeBuilder();  //

       HalfedgeDS he=builder->createMesh(V.rows(), F);  //

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
    
//    for(int i=0; i<neighbTri.size(); i++)
//    {   for(int j=0; j<neighbTri.at(i).size();j++)
//    { adjTris(i,j)=neighbTri.at(i).at(j);
////        cout<<neighbTri.at(i).at(j)<<" ";
////        cout<<"\n";
//    }
//    }
//    cout<<"size: "<< adjTris<<"\n";
//    vector<vector<double> > A;
    //   adjacency_list(F,A);
    
    
    
    
    regTri = MatrixXi::Zero(F.rows(),2); //1st col for marking&no of entries 2nd col for cluster id //for checking proxy assignment
    Vector3d clusterXi[k];
    Vector3d clusterNi[k];
    int arr[6] ={157, 349,23,308,130,422};
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
            clusterXi[i] = computeXi(ran);
            clusterNi[i] = computeNi(ran);
            
            for(int j=0; j<neighbTri.at(ran).size(); j++)
            {
//                cout<<"cent "<<clusterXi[i].transpose()<<"\t"<<"nor "<< clusterNi[i].transpose()<<"\n";
//                cout<<"before "<<i<<" "<<ran<<" "<<neighbTri.at(ran).at(j)<<endl;
                double l=l2error( clusterXi[i],clusterNi[i],neighbTri.at(ran).at(j));
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
                double l=l2error( clusterXi[tr.cluster],clusterNi[tr.cluster],neighbTri.at(tr.face).at(j));
                
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
//  std::cout<<
//    "Press '1' for per-face normals calling pre-defined functions of LibiGL."<<std::endl<<
//	"Press '2' for per-vertex normals calling pre-defined functions of LibiGL."<<std::endl<<
//    "Press '3' for lib_per-vertex normals using face-vertex structure of LibiGL ."<<std::endl<<
//	"Press '4' for HE_per-vertex normals using HalfEdge structure."<<std::endl;

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

//  Z colors
   // Use the z coordinate as a scalar field over the surface
    Z.col(0)=V.col(1);
//    for(int i=0; i<V.rows(); i++)
//    {
//        Z(i,0)=regTri(he.getFace(he.getEdge(i)),1);
//    }
//    for(int i=0; i<F.rows(); i++)
//    {
//        Z(i,0)=regTri(i,1);
//    }
    MatrixXd C;
  // Assign per-vertex colors
	igl::jet(Z,true,C);
	viewer.data().set_colors(C);  // Add per-vertex colors
  //viewer.core(0).align_camera_center(V, F);  //not needed
//    viewer.append_mesh();
	viewer.launch(); // run the viewer
    
//    clusters = bestfittingRegions( viewer, clusters, he , threshold, regTri);
}
