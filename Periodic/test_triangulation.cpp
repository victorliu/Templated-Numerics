#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath> // for ceil

#define PERIODIC_TRIANGULATION_DEBUG
#include "Triangulation.hpp"


typedef Periodic::Triangulation2<double> Triangulation;
typedef Triangulation::Lattice Lattice;
typedef Triangulation::UVCoord UVCoord;
typedef Triangulation::Pt2 Pt2;
typedef Triangulation::Offset Offset;
typedef Triangulation::Vec2 Vec2;
typedef Triangulation::LatticePoint LatticePoint;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Edge Edge;
typedef Triangulation::Face Face;
typedef Triangulation::ParticularEdge ParticularEdge;
typedef Triangulation::ParticularFace ParticularFace;


#define TPOSTSCRIPT_LIGHT
#include <AnalyticGeometry/TPostScript.hpp>

std::ostream& operator<<(std::ostream &os, const Triangulation &T){
	Vec2 u,v;
	T.GetLattice().GetBasis(u,v);
	TPostScript::DrawArrow(Pt2::Origin, u, os);
	TPostScript::DrawArrow(Pt2::Origin, v, os);

	std::stringstream buf;
	for(size_t i = 0; i < T.Vertices().size(); ++i){
		TPostScript::SetColor(0,0,0, os);
		TPostScript::DrawPoint(T.GetPoint(i), os);
		
		TPostScript::SetColor(1,0,0, os);
		buf.str("");
		buf << i;
		TPostScript::DrawText(T.GetPoint(i), buf.str(), os);
	}
	for(size_t i = 0; i < T.Faces().size(); ++i){
		size_t v[3];
		for(size_t j = 0; j < 3; ++j){
			v[j] = T.Faces()[i].v[j];
		}
		Pt2 p[3];
		T.GetFacePoints(i, p);
		TPostScript::SetColor(0,0,0, os);
		for(size_t j = 0; j < 3; ++j){
			//if(v[j] < v[(j+1)%3]){
				TPostScript::DrawSegment(p[j], p[(j+1)%3], os);
			//}
		}
		
		TPostScript::SetColor(0,0,1, os);
		buf.str("");
		buf << i;
		TPostScript::DrawText(p[0] + ((p[1]-p[0])+(p[2]-p[0]))/3., buf.str(), os);
	}
	return os;
}

void DumpTriangulation(std::ostream &os, const Triangulation &T){
	for(size_t i = 0; i < T.Vertices().size(); ++i){
		os << "V[" << i << "] = ("
			<< T.Vertices()[i].p.st[0] << "," << T.Vertices()[i].p.st[1]
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
}

int main(){
	std::ofstream f;
	
	Lattice L(Lattice::Vec2(1,0), Lattice::Vec2(0,1));
	Triangulation::vertex_vector_t vertices;
	vertices.push_back(Vertex(UVCoord(0.3,0.4)));
	
	vertices.push_back(Vertex(UVCoord(0.1,0.3)));
	vertices.push_back(Vertex(UVCoord(0.4,0)));
	vertices.push_back(Vertex(UVCoord(0.6,0.5)));
	
	Triangulation T(L, vertices);
	/*
	LatticePoint lp(UVCoord(0.1,0.1), Offset(-1,-1));
	ParticularFace face_info = T.ContainingFace(lp);
	std::cout << "Point at st(" << lp.st[0] << "," << lp.st[1] << ") is located in face " << face_info.face << "(" << face_info.offset[0] << "," << face_info.offset[1] << ")" << std::endl;
	*/
	
	/*
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
	*/
	f.open("out1.ps");
	TPostScript::Initialize<float>(f);
	f << T;
	TPostScript::Close<float>(f);
	f.close();
	
	//T.AddVertex(Vertex(Pt2(0.75,0.25)));
	//T.Remove(3);
	
	f.open("out2.ps");
	DumpTriangulation(f, T);
	f.close();
	
	return 0;
}
