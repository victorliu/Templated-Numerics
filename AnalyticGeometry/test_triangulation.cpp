#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TTriangulation2.hpp"

#define OUTPUT_POSTSCRIPT

typedef TTriangulation2<double> Triangulation;
typedef Triangulation::Pt2 Pt2;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Edge Edge;
typedef Triangulation::Face Face;


#ifdef OUTPUT_POSTSCRIPT

#define TPOSTSCRIPT_LIGHT
#include "TPostScript.hpp"

std::ostream& operator<<(std::ostream &os, const Triangulation &T){
	for(size_t i = 0; i < T.Vertices().size(); ++i){
		TPostScript::DrawPoint(T.Vertices()[i].p, os);
	}
	for(size_t i = 0; i < T.Faces().size(); ++i){
		if(!T.Faces()[i].IsInfinite()){
			for(size_t j = 0; j < 3; ++j){
				size_t v0 = T.Faces()[i].v[j];
				size_t v1 = T.Faces()[i].v[(j+1)%3];
				if(v0 < v1){
					TPostScript::DrawSegment(T.Vertices()[v0].p, T.Vertices()[v1].p, os);
				}
			}
		}else{
			size_t v0 = T.Faces()[i].v[0];
			size_t v1 = T.Faces()[i].v[1];
			TPostScript::DrawSegment(T.Vertices()[v0].p, T.Vertices()[v1].p, os);
		}
	}
	return os;
}

#else // OUTPUT_POSTSCRIPT

std::ostream& operator<<(std::ostream &os, const Triangulation &T){
	for(size_t i = 0; i < T.Vertices().size(); ++i){
		os << "V[" << i << "] = ("
			<< T.Vertices()[i].p.r[0] << "," << T.Vertices()[i].p.r[1]
			<< ")f("
			<< T.Vertices()[i].face
			<< ")" << std::endl;
	}
	for(size_t i = 0; i < T.Faces().size(); ++i){
		os << "F[" << i << "] = ("
			<< T.Faces()[i].v[0] << "," << T.Faces()[i].v[1] << "," << (int)(T.Faces()[i].v[2])
			<< ")n("
			<< T.Faces()[i].n[0] << "," << T.Faces()[i].n[1] << "," << T.Faces()[i].n[2]
			<< ")" << std::endl;
	}
	return os;
}

#endif // OUTPUT_POSTSCRIPT

int main(){
	std::ofstream f;
	
	Triangulation::vertex_vector_t vertices;
	vertices.push_back(Pt2(0,0));
	vertices.push_back(Pt2(0,1));
	vertices.push_back(Pt2(1,0));
	vertices.push_back(Pt2(1,1));
	Triangulation T(vertices);
	
	srand(0);
	
	for(size_t i = 0; i < 30; ++i){
		int p = rand()%100;
		if(p < 10){
			// Can't do this; not Delaunay-preserving
			// T.Split(Edge(rand() % T.Faces().size(), rand()%3));
		}else if(p < 25){
			if(T.Vertices().size() > 3){
				T.Remove(rand() % T.Vertices().size());
			}
		}else{
			Triangulation::Vertex v(Pt2(double(rand())/RAND_MAX, double(rand())/RAND_MAX));
			T.AddVertex(v);
		}
	}
	
	f.open("out1.ps");
	f << T;
	f.close();
	
	//T.AddVertex(Vertex(Pt2(0.75,0.25)));
	//T.Remove(3);
	/*
	f.open("out2.ps");
	f << T;
	f.close();
	*/
	return 0;
}
