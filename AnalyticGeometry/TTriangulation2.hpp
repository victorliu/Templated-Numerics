#ifndef _TTRIANGULATION2_HPP_
#define _TTRIANGULATION2_HPP_

#include "TPt2.h"
#include "TPredicates2.hpp"
#include <vector>
#include <list>

// Delaunay triangulation in 2D.

// Internally, only combinatoric structures are kept.
// There is a fictitious vertex at infinite, and all convex hull vertices
// are connected to it. The triangulation is therefore topologically
// equivalent to a sphere.

// Preprocessor switches:
//   TTRIANGULATION2_DEBUG - turns on some debug messages to std::cout

template <typename NumericType>
class TTriangulation2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	static const size_t infinity = size_t(-1);
	static const size_t invalid_face = size_t(-1);
	
	// The basic primtive elements of the triangulation: Vertex, Face and Edge
	struct Vertex{
		Pt2 p;       // actual location of vertex
		int tag;     // this is just some tag you can use for ID purposes
		// Don't mess with this; this is used by the triangulation to keep track of adjacency.
		size_t face;
		
		// We don't set face; face will be set when we do stuff to the vertices within this class
		Vertex(const Pt2 &P):p(P),tag(0){}
		Vertex(const Pt2 &P, int Tag):p(P),tag(Tag){}
	};
	struct Face{ // if one vertex is infinity, it is always v[2]
		size_t v[3]; // incident vertices, in positive orientation
		size_t n[3]; // neighboring faces across from corresponding vertex in v
		Face(size_t v0, size_t v1, size_t v2, size_t n0, size_t n1, size_t n2){
			v[0] = v0; v[1] = v1; v[2] = v2;
			n[0] = n0; n[1] = n1; n[2] = n2;
		}
		Face(const Face &f){
			v[0] = f.v[0]; v[1] = f.v[1]; v[2] = f.v[2];
			n[0] = f.n[0]; n[1] = f.n[1]; n[2] = f.n[2];
		}
		bool IsInfinite() const{ return (infinity == v[2]); }
		bool FaceAcrossFrom(size_t vertex_index, size_t &face_index) const{
			for(size_t i = 0; i < 3; ++i){
				if(v[i] == vertex_index){
					face_index = n[i];
					return true;
				}
			}
			return false;
		}
		int VertexOffset(size_t vertex_index) const{ // returns -1 if not found
			for(int i = 0; i < 3; ++i){
				if(vertex_index == v[i]){ return i; }
			}
			return -1;
		}
	};
	struct Edge{ // A combinatorial edge (each edge can be represented in two ways)
		size_t face; // an incident face
		int across;  // the vertex offset in Face.v of the face across which the edge sits (0-2)
		
		Edge(size_t incident_face, int opposite_vertex):face(incident_face),across(opposite_vertex%3){}
		Edge NextEdge() const{ return Edge(face, across+1); } // next edge of the face this belongs to (in the CCW direction)
		Edge PrevEdge() const{ return Edge(face, across+2); }
	};
	
	typedef std::vector<Vertex> vertex_vector_t;
	typedef std::vector<Face> face_vector_t;
	
public: // member functions
	// At all times we must have at least 3 vertices in the triangulation
	// If existing_vertices contains fewer than 3 vertices, the resulting
	// triangulation will be invalid. There is no way of recovering from this.
	TTriangulation2(std::vector<Vertex> &existing_vertices);
	TTriangulation2(const TTriangulation2& t);

	// Get sizes
	size_t NumVertices() const{ return vertices.size(); }
	size_t NumFaces() const{ return faces.size(); } // includes infinite faces
	size_t NumEdges() const{ return 3*NumVertices()/2; }
	size_t NumTriangles() const; // returns the number of finite faces
	
	// Get arrays
	const std::vector<Vertex>& Vertices() const{ return vertices; }
	const std::vector<Face>& Faces() const{ return faces; }
	
	// Find the face in which point p is located.
	// This always returns a face since we include infinite faces.
	size_t ContainingFace(const Pt2 &p, size_t face_index = 0);
	size_t AddVertex(const Vertex &v);
	// Give this a range of vertices
	template <class InputIterator>
	size_t AddVertices(InputIterator begin, InputIterator end);
	
	// Flip the edge across to the other diagonal
	void Flip(const Edge &e);
	
	// Note that this does not guarantee Delaunay-ness;
	// no other Delaunay functions may be used after calling this.
	size_t Split(const Edge &e, int tag = 0);
	
	// Remove a vertex and re-triangulate the hole
	bool Remove(size_t vertex_index);
	
	// End of public functions

private: // member data
	vertex_vector_t vertices;
	face_vector_t faces;
	typedef std::list<Edge> Hole;
private: // private helper functions

	//// Combinatorial primitives
	size_t AddVertexInFace(const Vertex &nv, size_t face_index){
		if(faces[face_index].IsInfinite()){
			return AddVertexInInfiniteFace(nv, face_index);
		}else{
			return AddVertexInInteriorFace(nv, face_index);
		}
	}
	size_t AddVertexInInteriorFace(const Vertex &nv, size_t face_index){
		// We assume that v lies in faces[face_index]
		
		// break up the face into 3 faces, and flip until good
		// original face's offset 0 vertex becomes the new vertex
		// Note that if v2 was infinity, all new infinite faces (just f1) have infinity at v2
		
		size_t v = vertices.size();
		vertices.push_back(nv); vertices.back().face = face_index;
		
		// The original face's vertices
		size_t v0 = faces[face_index].v[0];
		size_t v1 = faces[face_index].v[1];
		size_t v2 = faces[face_index].v[2];
		
		// These are the edges adjacent to v0 (v0-v1) and (v2-v0)
		Edge e1(EquivalentEdge(Edge(face_index, 1)));
		Edge e2(EquivalentEdge(Edge(face_index, 2)));
		
		// Make the new faces (these are completely set here)
		size_t f1 = faces.size();
		faces.push_back(Face(v0,v,v2, face_index,e1.face,infinity));
		size_t f2 = faces.size();
		faces.push_back(Face(v0,v1,v, face_index,infinity,e2.face));
		faces[f1].n[2] = f2;
		faces[f2].n[1] = f1;

		// Update neighbors of neighbors of face_index
		faces[e1.face].n[e1.across] = f1;
		faces[e2.face].n[e2.across] = f2;

		// Finally fix up this face
		faces[face_index].v[0] = v;
		faces[face_index].n[1] = f1;
		faces[face_index].n[2] = f2;
		
		if(vertices[v0].face == face_index){
			vertices[v0].face = f2;
		}

		return v;
	}
	Edge EquivalentEdge(const Edge &e) const{ // given an edge, return an equivalent representation of the edge using the face opposite the edge
		size_t other_face = faces[e.face].n[e.across];
		int offset = faces[other_face].VertexOffset(faces[e.face].v[(e.across+1)%3]);
		return Edge(other_face, (offset+1)%3);
	}
	void MakeHole(size_t vertex_in_hole, Hole &hole){
#ifdef TTRIANGULATION2_DEBUG
		std::cout << "In MakeHole, removing faces around v" << vertex_in_hole << std::endl;
#endif
		// Note that this function does not delete the vertex; it only removes all faces around a vertex
		std::vector<size_t> faces_in_hole;
		
		size_t cur_face = vertices[vertex_in_hole].face;
		size_t start = cur_face;
		Edge e(cur_face, (faces[cur_face].VertexOffset(vertex_in_hole)+2)%3);
		do{
			//std::cout << "Iteration edge between vertices " << faces[e.face].v[(e.across+1)%3] << "," << faces[e.face].v[(e.across+2)%3] << std::endl;
			// Move the hole boundary vertices' face away from the faces that will be deleted
			Edge en(EquivalentEdge(e.NextEdge()));
			//std::cout << "Neighbor edge between vertices " << faces[en.face].v[(en.across+1)%3] << "," << faces[en.face].v[(en.across+2)%3] << std::endl;
			size_t v;
			v = faces[cur_face].v[(en.across+1)%3];
			if(infinity != v && vertices[v].face == cur_face){ vertices[v].face = en.face; }
			v = faces[cur_face].v[(en.across+2)%3];
			if(infinity != v && vertices[v].face == cur_face){ vertices[v].face = en.face; }
			
			// Change the neighbors to indicate a deletion
			faces[en.face].n[en.across] = invalid_face;
			
			hole.push_back(en);
			faces_in_hole.push_back(cur_face);
			
			e = EquivalentEdge(e.PrevEdge());
			cur_face = e.face;
		}while(cur_face != start);
		
		std::sort(faces_in_hole.begin(), faces_in_hole.end());

#ifdef TTRIANGULATION2_DEBUG
		for(size_t i = 0; i < faces_in_hole.size(); ++i){
			std::cout << faces_in_hole[i] << " ";
		}std::cout << std::endl;
#endif

		// Update all indices
		// Find the mapping from old face indices to new indices
		std::vector<size_t> new_face_indices(faces.size());
		size_t j = 0;
		for(size_t i = 0; i < faces.size(); ++i){
			if(i > faces_in_hole[j]){
				++j;
			}
			new_face_indices[i] = i-j;
		}
		// Remove all the faces
		for(size_t i = 0; i < faces_in_hole.size(); ++i){
			faces.erase(faces.begin()+faces_in_hole[i]-i);
		}
		// Update all the places where there are face indices
		for(size_t i = 0; i < vertices.size(); ++i){
			if(invalid_face != vertices[i].face){
				vertices[i].face = new_face_indices[vertices[i].face];
			}
		}
		for(size_t i = 0; i < faces.size(); ++i){
			for(size_t j = 0; j < 3; ++j){
				if(invalid_face != faces[i].n[j]){
					faces[i].n[j] = new_face_indices[faces[i].n[j]];
				}
			}
		}
		for(typename Hole::iterator i = hole.begin(); i != hole.end(); ++i){
			if(invalid_face != i->face){
				i->face = new_face_indices[i->face];
			}
		}
	}
	
	//// Specifically Delaunay primitives
	void FillHole(const Hole &hole_to_fill){
#ifdef TTRIANGULATION2_DEBUG
		std::cout << "In FillHole, hole:" << std::endl;
		for(typename Hole::const_iterator i = hole_to_fill.begin(); i != hole_to_fill.end(); ++i){
			std::cout << "v" << (int)faces[i->face].v[(i->across+1)%3] << "-v" << (int)faces[i->face].v[(i->across+2)%3] << " ";
		}std::cout << std::endl;
#endif		
		// Fills in a hole using a Delaunay triangulation
		// The hole is described by a list of bounding edges,
		// in positive orientation (going around the hole in CCW order)
		
		if(hole_to_fill.size() < 3){ return; } // invalid hole
		
		// We know ahead of time how many new faces we will add
		faces.reserve(faces.size() + hole_to_fill.size() - 2);
		
		std::list<Hole> holes;
		holes.push_back(hole_to_fill);
		
		while(holes.size() > 0){
			Hole hole = holes.front();
			holes.pop_front();
			
			if(3 == hole.size()){
				// Create a face and we're done
				size_t new_face_index = faces.size();
				typename Hole::const_iterator h = hole.begin();
				size_t v[3], n[3];
				int inf_offset = 0;
				for(size_t i = 0; i < 3; ++i){
					n[i] = h->face;
					v[i] = faces[h->face].v[(h->across+1)%3];
					if(infinity == v[i]){ inf_offset = i; }
					faces[h->face].n[h->across] = new_face_index; // update neighbors preemptively
					++h;
				}
				faces.push_back(Face(v[(inf_offset+1)%3],v[(inf_offset+2)%3],v[inf_offset], n[inf_offset],n[(inf_offset+1)%3],n[(inf_offset+2)%3]));
				continue;
			}
			
			// If the hole has infinite boundary edges, make sure they're not at the beginning
			// We cycle through the hole's edges until the first one is not infinite
			while(true){
				size_t f = hole.front().face;
				int a = hole.front().across;
				if((infinity == faces[f].v[(a+1)%3]) || (infinity == faces[f].v[(a+2)%3])){
					hole.push_back(hole.front());
					hole.pop_front();
				}else{ break; }
			}
			// Note that the above cyclic permutation always terminates, as there is only
			// one infinite vertex, and the size of the hole is at least 3.
			
			// We will now add a face with one of its edges being e
			Edge e = hole.front(); hole.pop_front();
			size_t v0 = faces[e.face].v[(e.across+2)%3];
			size_t v1 = faces[e.face].v[(e.across+1)%3];
			size_t v2 = infinity;
			// The new triangle will be v0,v1,v2, where v2 is to be determined
			size_t v_last_non_inf; // This is for recording state from the following iteration
			const Pt2& p0 = vertices[v0].p;
			const Pt2& p1 = vertices[v1].p;
			
			// Note that if the new face must have one vertex at infinity, it will be correctly place in v2!
			
			// For all the other hole edges, we must now find the third vertex to join with it
			// We iterate over all the possible third vertices
			typename Hole::iterator hend = hole.end(); --hend; // the last vertex is on the edge one before the edge before e
			typename Hole::iterator h = hole.begin(); // hole.begin() is actually the edge after e since we popped away e
			typename Hole::iterator cut_location(h); // we will be cutting hole after this iterator location
			while(h != hend){
				size_t cur_vertex = faces[h->face].v[(h->across+1)%3];
				// cur_vertex will take on all the proper possible values for v2 along the hole boundary
				// We will go through all of them without early exit to get the last one
				
				// Throughout this iteration, v2 will store the last (in
				// the sense of previous) possible vertex to use.
				// Also, v_last_non_inf is used to store the last
				// possible non-infinite vertex that we can use.
				
				// Now we only have to check the Delaunay criterion.
				if(infinity == cur_vertex){
					if(infinity == v2){ // If we have not yet found a finite vertex, then consider putting v0,v1 on the hull
						cut_location = h;
					} // If we HAVE found a finite vertex, we cannot form a triangle towards infinity since we would not be forming a convex hull.
				}else{ // v0 or v1 are not infinite due to the cyclic permutation above
					const Pt2 &p = vertices[cur_vertex].p;
					if(Orient2(p0, p1, p) > 0){ // if the proposed face has the proper orientation...
						if(
							(infinity == v2) || // If we have not yet found a finite vertex, update v2 (the last possible vertex) to this finite one
							(InCircle2(p0, p1, vertices[v_last_non_inf].p, p) > 0) // If p is INSIDE the last circumcircle, then we better use p instead!
						){
							v2 = cur_vertex;
							v_last_non_inf = cur_vertex;
							cut_location = h;
						}
					}
				}
				
				++h;
			}
			
			// We now need to make a face defined by edge e and vertex v2
			// We also need to split the hole into two holes since the triangle could have
			// not clipped an ear off the hole.
			
			// First we check for the two cases where the hole is not actually split.
			Edge en(e); // the neighboring edge directly after/before e (we have to initialize it with something, so use e)
			int offset;
			if( // assign en to edge directly after e (in CCW order) and see if the far vertex is v2
				((en = hole.front()), (v2 == faces[en.face].v[(en.across+1)%3]))
			){ // If v2 is directly after v0 and v1, then there is no hole split
#ifdef TTRIANGULATION2_DEBUG
				std::cout << " In case 1; no split" << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v0,v1,v2, en.face,invalid_face,e.face));
				// Update neighbors
				faces[e.face].n[e.across] = new_face_index;
				faces[en.face].n[en.across] = new_face_index;
				
				// The hole's old edge needs replacing with the new face's edge
				hole.pop_front();
				hole.push_front(Edge(new_face_index, 1));
				holes.push_front(hole);
			}else if( // assign en to edge directly before e (in CCW order) and see if the far vertex is v2
				((en = hole.front()), (v2 == faces[en.face].v[(en.across+2)%3]))
			){ // If v2 is directly before v0 and v1, then there is no hole split
#ifdef TTRIANGULATION2_DEBUG
				std::cout << " In case 2; no split" << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v0,v1,v2, invalid_face,en.face,e.face));
				// Update neighbors
				faces[e.face].n[e.across] = new_face_index;
				faces[en.face].n[en.across] = new_face_index;
				
				// The hole's old edge needs replacing with the new face's edge
				hole.pop_front();
				hole.push_front(Edge(new_face_index, 0));
				holes.push_front(hole);
			}else{ // We will have a hole split
#ifdef TTRIANGULATION2_DEBUG
				std::cout << " Splitting hole" << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v0,v1,v2, invalid_face,invalid_face,e.face));
				faces[e.face].n[e.across] = new_face_index; // update neighbor
				
				// Split the hole
				Hole new_hole;
				++cut_location;
				while(hole.begin() != cut_location){
					new_hole.push_back(hole.front()); // note this has to be push_back to keep the orientation the same
					hole.pop_front();
				}
				hole.push_front(Edge(new_face_index,1));
				new_hole.push_front(Edge(new_face_index,0));
				holes.push_front(hole);
				holes.push_front(new_hole);
			}
		}
	}
	
	size_t AddVertexInInfiniteFace(const Vertex &nv, size_t face_index){
		// Assumes nv is located in faces[face_index]
		// infinite vertex is always at faces[face_index].v[2]
		
		// Find all neighboring infinite faces that will need to be modified
		// Only a set of contiguous neighbors need modification since the
		// hull is convex and they are the only visible faces from nv.
		std::vector<size_t> inf_faces_before, inf_faces_after;
		// Find those faces before (going clockwise around the infinite point)
		{
			Edge e(EquivalentEdge(Edge(face_index,1)));
			while(true){
				if(Orient2(nv.p, vertices[faces[e.face].v[0]].p, vertices[faces[e.face].v[1]].p) > 0){
					inf_faces_before.push_back(e.face);
				}else{ break; }
				e.across = (e.across+2)%3;
				e = EquivalentEdge(e);
			}
		}
		// Find those faces after (going counter-clockwise around the infinite point)
		{
			Edge e(EquivalentEdge(Edge(face_index,0)));
			while(true){
				if(Orient2(nv.p, vertices[faces[e.face].v[0]].p, vertices[faces[e.face].v[1]].p) > 0){
					inf_faces_after.push_back(e.face);
				}else{ break; }
				e.across = (e.across+1)%3;
				e = EquivalentEdge(e);
			}
		}
		
		// Now we can actually add the vertex
		size_t v = AddVertexInInteriorFace(nv, face_index);
		
		// Add all the other faces (these form very coney shapes)
		for(size_t i = 0; i < inf_faces_before.size(); ++i){
			Flip(Edge(inf_faces_before[i],0));
		}
		for(size_t i = 0; i < inf_faces_after.size(); ++i){
			Flip(Edge(inf_faces_after[i],1));
		}
		
		return v;
	}
	// Make this vertex locally Delaunay; keep flipping if necessary
	void MakeDelaunay(size_t vertex_index){
		size_t f = vertices[vertex_index].face;
		size_t next, start = f;
		do{
			int i = faces[f].VertexOffset(vertex_index);
			next = faces[f].n[(i+1)%3];
			PropagatingDelaunayFlip(Edge(f, i));
			f = next;
		}while(next != start);
	}
	// Check if an edge needs flipping, if so, keep propagating the flips
	void PropagatingDelaunayFlip(const Edge &e){
		Edge e2(EquivalentEdge(e));
		if(InCircumcircle(e2.face, vertices[faces[e.face].v[e.across]].p) < 0){
			return;
		}
		Flip(e);
		// Possibly flip the two edges of the neighboring face adjacent to edge e
		PropagatingDelaunayFlip(e);
		PropagatingDelaunayFlip(Edge(e2.face, faces[e2.face].VertexOffset(faces[e.face].v[e.across])));
	}
	NumericType InCircumcircle(size_t face_index, const Pt2 &p){
		if(!faces[face_index].IsInfinite()){
			return InCircle2(
				vertices[faces[face_index].v[0]].p,
				vertices[faces[face_index].v[1]].p,
				vertices[faces[face_index].v[2]].p,
				p);
		}else{
			return Orient2(
				vertices[faces[face_index].v[0]].p,
				vertices[faces[face_index].v[1]].p,
				p);
		}
	}
private:
	// Embedded GameRand implementation for location queries
	unsigned int rng_low;
	unsigned int rng_high;
	unsigned int rng(){
		rng_high = (rng_high << 16) + (rng_high >> 16);
		rng_high += rng_low; rng_low += rng_high;
		return rng_high;
	}
};


// At all times we must have at least 3 vertices in the triangulation
template <typename NumericType>
TTriangulation2<NumericType>::TTriangulation2(std::vector<Vertex> &existing_vertices){
	// initialize GameRand
	rng_high = 0xDEADBEEF;
	rng_low = rng_high ^ 0x49616E42;
	
	if(existing_vertices.size() >= 3){
		vertices.reserve(existing_vertices.size());
		// Insert first face
		vertices.insert(vertices.end(), existing_vertices.begin(), existing_vertices.begin()+3);
		if(Orient2(vertices[0].p, vertices[1].p, vertices[2].p) > 0){
			faces.push_back(Face(0,1,2, 2,3,1));
			faces.push_back(Face(1,0,infinity, 3,2,0));
			faces.push_back(Face(2,1,infinity, 1,3,0));
			faces.push_back(Face(0,2,infinity, 2,1,0));
		}else{
			faces.push_back(Face(0,2,1, 2,3,1));
			faces.push_back(Face(2,0,infinity, 3,2,0));
			faces.push_back(Face(1,2,infinity, 1,3,0));
			faces.push_back(Face(0,1,infinity, 2,1,0));
		}
		vertices[0].face = 0;
		vertices[1].face = 0;
		vertices[2].face = 0;
		// Insert remaining faces
		AddVertices(existing_vertices.begin()+3,existing_vertices.end());
	}
}

template <typename NumericType>
TTriangulation2<NumericType>::TTriangulation2(const TTriangulation2& t):
	vertices(t.vertices),
	faces(t.faces),
	rng_low(t.rng_low),
	rng_high(t.rng_high)
{}

template <typename NumericType>
size_t TTriangulation2<NumericType>::NumTriangles() const{
	size_t total = 0;
	for(typename face_vector_t::const_iterator i = faces.begin(); i != faces.end(); ++i){
		if(!(i->IsInfinite())){
			++total;
		}
	}
	return total;
}

template <typename NumericType>
size_t TTriangulation2<NumericType>::ContainingFace(const Pt2 &p, size_t face_index = 0){
	if(faces[face_index].IsInfinite()){
		faces[face_index].FaceAcrossFrom(infinity, face_index);
	}
	while(true){
		if(faces[face_index].IsInfinite()){
			return face_index;
		}
		int outside[2]; int noutside = 0;
		for(size_t i = 0; i < 3; ++i){
			if(Orient2(vertices[faces[face_index].v[i]].p, vertices[faces[face_index].v[(i+1)%3]].p, p) < 0){
				outside[noutside++] = faces[face_index].n[(i+2)%3];
				if(noutside >= 2){ break; }
			}
		}
		if(0 == noutside){
			return face_index;
		}else if(1 == noutside){
			face_index = outside[0];
		}else{
			// pick a random side to flip over to
			face_index = outside[rng()&1];
		}
	}
}

template <typename NumericType>
size_t TTriangulation2<NumericType>::AddVertex(const Vertex &v){
	size_t face_index = 0;
	face_index = ContainingFace(v.p, face_index);
	size_t ret = AddVertexInFace(v, face_index);
	MakeDelaunay(ret);
	return ret;
}
	
template <typename NumericType>
template <class InputIterator>
size_t TTriangulation2<NumericType>::AddVertices(InputIterator begin, InputIterator end){
	size_t face_index = 0;
	size_t ret = infinity;
	for(InputIterator i = begin; i != end; ++i){
		face_index = ContainingFace(i->p, face_index);
		size_t idx = AddVertexInFace(*i, face_index);
		MakeDelaunay(idx);
		if(ret != infinity){ ret = idx; }
	}
	return ret;
}

template <typename NumericType>
void TTriangulation2<NumericType>::Flip(const Edge &e){
	// faces[e.face].v[e.across] is maintained throughout this process
	
	Edge e2(EquivalentEdge(e)); // this is the same edge from the other face
	
	// This is the layout before we do anything
	//       e.face[e.across]
	//        |/ e_next |/
	//   -----+---------+ -- v_prev
	//       /|        /|
	//      / | e.face/ |
	//     /  |      /  |
	//        |     /   |
	//        |   e/e2  |
	//        |   /     |
	//        |  /      |
	//        | /e2.face| /
	//        |/        |/
	//  ------+---------+ -- faces[n].v[nacross]
	//        | e2_next/|
	//       v_next
	size_t v_prev = faces[e.face].v[(e.across+2)%3];
	size_t v_next = faces[e.face].v[(e.across+1)%3];

	Edge  e_next(EquivalentEdge( e.NextEdge()));
	Edge e2_next(EquivalentEdge(e2.NextEdge()));

	// The faces get the edge flip
	faces[ e.face].v[( e.across+2)%3] = faces[e2.face].v[e2.across];
	faces[e2.face].v[(e2.across+2)%3] = faces[ e.face].v[ e.across];

	// Update the neighbors
	faces[e.face].n[e.across] = e2_next.face;
	faces[e2_next.face].n[e2_next.across] = e.face;

	faces[ e.face].n[( e.across+1)%3] = e2.face;
	faces[e2.face].n[(e2.across+1)%3] =  e.face;

	faces[e2.face].n[e2.across] = e_next.face;
	faces[e_next.face].n[e_next.across] = e2.face;

	// Update incident face information
	if(vertices[v_prev].face == e.face){
		vertices[v_prev].face = e2.face;
	}
	if(vertices[v_next].face == e2.face){
		vertices[v_next].face = e.face;
	}
}

template <typename NumericType>
size_t TTriangulation2<NumericType>::Split(const Edge &e, int tag = 0){ // returns new vertex index
	 // cannot split an infinite edge
	if(faces[e.face].IsInfinite()){
		if(2 != e.across){ // if e.across == 2, the edge is on convex hull
			return infinity;
		}
	}
	Edge e2(EquivalentEdge(e));
	Pt2 p0(vertices[faces[e.face].v[(e.across+1)%3]].p);
	Pt2 p1(vertices[faces[e.face].v[(e.across+1)%3]].p);
	p0 += (NumericType(1)/NumericType(2)) * (p1-p0);
	size_t v = AddVertexInFace(Vertex(p0, tag), e.face);
	Flip(e2);
	return v;
}

template <typename NumericType>
bool TTriangulation2<NumericType>::Remove(size_t vertex_index){
	if(vertices.size() < 4){ return false; }
	Hole hole;
	MakeHole(vertex_index, hole);
	FillHole(hole);
	
	for(size_t i = 0; i < faces.size(); ++i){
		for(size_t j = 0; j < 3; ++j){
			if(infinity != faces[i].v[j] && faces[i].v[j] > vertex_index){ --faces[i].v[j]; }
		}
	}
	vertices.erase(vertices.begin()+vertex_index);
	return true;
}


#endif // _TTRIANGULATION2_HPP_
