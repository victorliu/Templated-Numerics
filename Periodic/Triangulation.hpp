#ifndef _PERIODIC_TRIANGULATION_HPP_
#define _PERIODIC_TRIANGULATION_HPP_

// NONWORKING
// TODO: FIX AddVertex

#include <vector>
#include <list>
#include <stack>
#include <set>
#include "Lattice2.hpp"
#include <AnalyticGeometry/TTriangulation2.hpp>

// Periodic Delaunay triangulation in 2D.

// See AnalyticGeometry/Triangulation.hpp for details
// Throughout the code, UV coordinates are coefficients of the lattice
// basis vectors, so that the actual point is uv[0]*L.u+uv[1]*L.v
// We use the notation s-t coordinates interchangeably with UV coords.

// Preprocessor switches:
//   PERIODIC_TRIANGULATION_DEBUG - turns on some debug messages to std::cout

// This file proves again that there IS such a thing as too many comments.

namespace Periodic{

template <typename NumericType>
class Triangulation2{
public:
	typedef NumericType value_type;
	
	typedef Lattice2<NumericType> Lattice;
	typedef typename Lattice::Pt2 Pt2;
	typedef typename Lattice::Vec2 Vec2;
	typedef typename Lattice::Mat2 Mat2;
	typedef typename Lattice::Coord2 Coord2;
	typedef typename Lattice::UVCoord UVCoord;
	
	static const size_t invalid_face = size_t(-1);
	
	struct Offset{
		int off[2];
		Offset(){ off[0] = off[1] = 0; }
		Offset(int o1, int o2){ off[0] = o1; off[1] = o2; }
		Offset(int o[2]){ off[0] = o[0]; off[1] = o[1]; }
		Offset(const Offset &o){ off[0] = o.off[0]; off[1] = o.off[1]; }
		Offset& operator=(const Offset &o){ off[0] = o.off[0]; off[1] = o.off[1]; return *this; }
		Offset& operator+=(const Offset &o){ off[0] += o.off[0]; off[1] += o.off[1]; return *this; }
		Offset& operator-=(const Offset &o){ off[0] -= o.off[0]; off[1] -= o.off[1]; return *this; }
		bool operator==(const Offset &o) const{ return off[0] == o.off[0] && off[1] == o.off[1]; }
		int operator[](int i) const{ return off[i]; }
		int& operator[](int i){ return off[i]; }
		bool operator<(const Offset &o) const{ return ((off[0] < o.off[0]) || (off[0] == o.off[0] && off[1] < o.off[1])); }
	};
	struct LatticePoint{
		UVCoord st;    // Lattice uv coordinates of point, always in the range [0,1)^2
		Offset offset; // The offset of the fundamental parallelogram in which the vertex sits
		LatticePoint(const UVCoord &uv):st(uv),offset(0,0){}
		LatticePoint(const UVCoord &uv, const Offset &off):st(uv),offset(off){}
	};
	// The basic primtive elements of the triangulation: Vertex, Face and Edge
	struct Vertex{
		LatticePoint p; // the actual point where the vertex is located
		int tag;        // this is just some tag you can use for ID purposes
		// Don't mess with this; this is used by the triangulation to keep track of adjacency.
		size_t face;
		
		// We don't set face; face will be set when we do stuff to the vertices within this class
		Vertex(const UVCoord &uv):p(uv),tag(0){}
		Vertex(const UVCoord &uv, int Tag):p(uv),tag(Tag){}
		Vertex(const LatticePoint &P, int Tag):p(P),tag(Tag){}
		Vertex(const Vertex& v, const Offset &offset):p(v.p),tag(v.tag){
			p.offset = offset;
		}
	};
	struct Face{
		size_t v[3]; // incident vertices, in positive orientation
		size_t n[3]; // neighboring faces across from corresponding vertex in v
		
		// The vertices of a face may not lie in a single fundamental parallelogram (unit cell)
		// However, they may not lie in non-neighboring cells.
		// This means that each each vertex can only have a relative offset of +/-1 with
		// respect to any other vertex. Also, if the triangulation is not 1-sheeted,
		// all copies of the vertices with their corresponding offsets are maintained.
		// Therefore, a face need not store the actual offset of each vertex, only those
		// that "wrap".
		// For faces near the right or upper edges of the n-sheeted unit cell, each vertex
		// index obviously refers to a vertex within the n-sheeted unit cell, and the
		// corresponding offset value is 0 if the face does not wrap to the left/bottom
		// edge, or 1 if it does. In this way, we only store faces that wrap around
		// the right/top sides, and never those that wrap around the left/bottom sides
		// Only one bit is needed to encode the offset then, and we use a bit array:
		unsigned char wraps;
		// The bit packing of offsets is as follows:
		// 00tststs
		//      |||
		//      ||Offset of s coord of v[0]
		//      |Offset of t coord of v[0]
		//      Offset of s coord of v[1]
		//     etc.
		
		Face(size_t v0, size_t v1, size_t v2, size_t n0, size_t n1, size_t n2, unsigned char Wraps = 0):wraps(Wraps){
			v[0] = v0; v[1] = v1; v[2] = v2;
			n[0] = n0; n[1] = n1; n[2] = n2;
		}
		Face(const Face &f):wraps(f.wraps){
			v[0] = f.v[0]; v[1] = f.v[1]; v[2] = f.v[2];
			n[0] = f.n[0]; n[1] = f.n[1]; n[2] = f.n[2];
		}
		bool FaceAcrossFrom(size_t vertex_index, size_t &face_index) const{
			for(size_t i = 0; i < 3; ++i){
				if(v[i] == vertex_index){
					face_index = n[i];
					return true;
				}
			}
			return false;
		}
		// This returns the offset within the v array, not the periodic offset
		int VertexOffset(size_t vertex_index) const{ // returns -1 if not found
			for(int i = 0; i < 3; ++i){
				if(vertex_index == v[i]){ return i; }
			}
			return -1;
		}
		int GetVertexWrap(int i) const{
			i <<= 1;
			return ((wraps >> i) & 0x3);
		}
	};
	struct ParticularFace{
		size_t face;
		Offset offset; // offset of the face
		ParticularFace(const ParticularFace &f):face(f.face),offset(f.offset){}
		ParticularFace(size_t face_index, const Offset &off = Offset()):face(face_index),offset(off){}
		bool operator<(const ParticularFace &f) const{ return ((face < f.face) || (face == f.face && offset < f.offset)); }
	};
	struct Edge{ // A combinatorial edge (each edge can be represented in two ways)
		size_t face; // an incident face
		int across;  // the vertex offset in Face.v of the face across which the edge sits (0-2)
		
		Edge(size_t incident_face, int opposite_vertex):face(incident_face),across(opposite_vertex%3){}
		Edge NextEdge() const{ return Edge(face, across+1); } // next edge of the face this belongs to (in the CCW direction)
		Edge PrevEdge() const{ return Edge(face, across+2); }
		bool operator==(const Edge &e) const{ return face == e.face && across == e.across; }
	};
	struct ParticularEdge{
		Edge edge;
		Offset offset; // offset of the face
		ParticularEdge(const ParticularEdge &e):edge(e.edge),offset(e.offset){}
		ParticularEdge(const Edge &e, const Offset &off = Offset()):edge(e),offset(off){}
		ParticularEdge(size_t incident_face, int opposite_vertex, const Offset &off = Offset()):edge(incident_face,opposite_vertex),offset(off){}
		ParticularEdge(const ParticularFace &face_info, int opposite_vertex):edge(face_info.face,opposite_vertex),offset(face_info.offset){}
		ParticularFace GetParticularFace() const{ return ParticularFace(edge.face, offset); }
		ParticularEdge NextEdge() const{ return ParticularEdge(edge.NextEdge(), offset); } // next edge of the face this belongs to (in the CCW direction)
		ParticularEdge PrevEdge() const{ return ParticularEdge(edge.PrevEdge(), offset); }
		bool operator==(const ParticularEdge &e) const{ return edge == e.edge && offset == e.offset; }
	};
	// The Particular* structures refer to one of the copies within the covering domain
	
	typedef std::vector<Vertex> vertex_vector_t;
	typedef std::vector<Face> face_vector_t;
	
public: // member functions
	// At all times we must have at least 1 vertex in the triangulation
	// If existing_vertices is empty, the resulting triangulation will be invalid.
	// There is no way of recovering from this.
	// The offsets in existing_vertices are ignored.
	Triangulation2(const Lattice &lattice, std::vector<Vertex> &existing_vertices);
	Triangulation2(const Triangulation2& t);

	const Lattice& GetLattice() const{ return L; }

	// Get sizes
	size_t NumVertices() const{ return vertices.size(); }
	size_t NumFaces() const{ return faces.size(); } // includes infinite faces
	size_t NumEdges() const{ return 3*NumVertices()/2; }
	size_t NumTriangles() const{ return faces.size(); }
	
	Pt2 GetPoint(size_t vertex_index) const;
	void GetFacePoints(const ParticularFace &face_info, Pt2 v[3]) const;
	Pt2 GetFacePoint(const ParticularFace &face_info, int which) const;
	
	// Get arrays
	const std::vector<Vertex>& Vertices() const{ return vertices; }
	const std::vector<Face>& Faces() const{ return faces; }
	
	// Find the face in which point p is located.
	// containing_face_offset is either zero or -(2*covering_extent+1) in each direction
	ParticularFace ContainingFace(const LatticePoint &p, size_t face_guess = 0) const;
	size_t AddVertex(const Vertex &v, size_t *face_index = NULL);
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
	Lattice L;
	vertex_vector_t vertices;
	face_vector_t faces;
	
	int covering_extent[2];
	bool one_sheeted;
	
	NumericType edge_length_threshold2; // square of threshold length; if all edges are shorter, we can be sure there are no self edges
	
	typedef std::map<size_t, std::list<size_t> > long_edge_map_t;
	long_edge_map_t long_edges; // map of vertex indices to neighboring vertex indices
	
	// vertex_to_copies contains as keys all zero-offset vertex indices.
	// It maps to a vector of all copies of the corresponding zero-offset vertex, including the zero offset one
	// Therefore, the length of the vector should be (2*covering_extent[0]+1)*(2*covering_extent[1]+1)
	typedef std::map<size_t, std::vector<size_t> > vertex_to_copies_map_t;
	vertex_to_copies_map_t vertex_to_copies;
	
	// copy_to_vertex contains as keys each vertex index, and maps to the indices
	// of the corresponding zero-offset vertex index.
	// All zero-offset vertices are also contained in this map, and they map to themselves.
	typedef std::map<size_t, size_t> copy_to_vertex_map_t;
	copy_to_vertex_map_t copy_to_vertex;
	
	typedef std::list<ParticularEdge> Hole;
private: // private helper functions

	//// Combinatorial primitives
	size_t AddVertexInFace(const Vertex &nv, ParticularFace face_info){ // we need a copy face_info since the face indices will change
		// We assume that v lies in faces[face_index]
		
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "In AddVertexInFace, new vertex at "
			<< nv.p.offset[0] << "+" << nv.p.st[0] << "," << nv.p.offset[1] << "+" << nv.p.st[1]
			<< " in face " << face_info.face
			<< " with offset " << face_info.offset[0] << "," << face_info.offset[1] << std::endl;
#endif

		// First, find all faces in conflict with the new vertex
		// It is assumed that faces[face_index] (with face_offset) is known to be in conflect
		// Recursively find all neighbors of the current face that are in conflict
		std::set<ParticularFace> faces_in_conflict;
		std::set<size_t> faces_in_conflict_indices;
		faces_in_conflict.insert(face_info);
		faces_in_conflict_indices.insert(face_info.face);
		Hole hole;
		{
			std::stack<ParticularEdge> edges_to_test;
			edges_to_test.push(EquivalentEdge(ParticularEdge(face_info, 0)));
			hole.push_back(edges_to_test.top());
			edges_to_test.push(EquivalentEdge(ParticularEdge(face_info, 1)));
			hole.push_back(edges_to_test.top());
			edges_to_test.push(EquivalentEdge(ParticularEdge(face_info, 2)));
			hole.push_back(edges_to_test.top());
			while(!edges_to_test.empty()){
				ParticularEdge cur_edge(edges_to_test.top()); edges_to_test.pop();
				ParticularFace cur_face(cur_edge.GetParticularFace());
				if(InCircumcircle(cur_face, nv.p) > 0){
#ifdef PERIODIC_TRIANGULATION_DEBUG
					std::cout << "   Hole: ";
					for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
						size_t v1 = faces[i->edge.face].v[(i->edge.across+2)%3];
						size_t v2 = faces[i->edge.face].v[(i->edge.across+1)%3];
						std::cout << "(" << v1 << "-" << v2 << ")o(" << i->offset[0] << "," << i->offset[1] << ") ";
					}std::cout << std::endl;
					std::cout << "   Conflicting face: " << cur_face.face << "(" << cur_face.offset[0] << "," << cur_face.offset[1] << ")" << std::endl;
#endif
					faces_in_conflict.insert(cur_face);
					faces_in_conflict_indices.insert(cur_face.face);
					// Expand the hole
					// Find cur_edge in the hole
					typename Hole::iterator i;
					for(i = hole.begin(); i != hole.end(); ++i){
						if((*i) == cur_edge){ break; }
					}
					// We should have i != hole.end()
					
					i = hole.erase(i); // after this, i points to the place after where we want to insert, perfect!
					
					size_t f1, f2; // store the next faces to check; for assigning vertex faces
					// Note the order of these two; must insert in reverse order to chain the iterator assignment
					edges_to_test.push(ParticularEdge(EquivalentEdge(cur_edge.PrevEdge())));
					i = hole.insert(i, edges_to_test.top()); f2 = edges_to_test.top().edge.face;
					edges_to_test.push(ParticularEdge(EquivalentEdge(cur_edge.NextEdge())));
					i = hole.insert(i, edges_to_test.top()); f1 = edges_to_test.top().edge.face;
				}
			}
		}
#ifdef PERIODIC_TRIANGULATION_DEBUG
		{
			std::cout << "   Faces in conflict: ";
			typename std::set<ParticularFace>::const_iterator i;
			for(i = faces_in_conflict.begin(); i != faces_in_conflict.end(); ++i){
				std::cout << i->face << "(" << i->offset[0] << "," << i->offset[1] << ") ";
			}
			std::cout << std::endl;
		}
#endif

		// Move the face of the vertices on the boundary of the hole away
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Updating vertices[i].face for i on boundary" << std::endl;
#endif		
		for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
			for(size_t j = 0; j < 2; ++j){
				size_t vj = faces[i->edge.face].v[(i->edge.across+j+1)%3];
				size_t old_face = vertices[vj].face;
				if(faces_in_conflict_indices.find(old_face) != faces_in_conflict_indices.end()){
#ifdef PERIODIC_TRIANGULATION_DEBUG
					std::cout << "      Looking for new face of vertex " << vj << ", current face: " << old_face << std::endl;
#endif
					// Circle around the vertex to find a new face
					Edge ve(old_face, faces[old_face].VertexOffset(vj)+1);
					do{
						ve = EquivalentEdge(ve).PrevEdge();
#ifdef PERIODIC_TRIANGULATION_DEBUG
						std::cout << "         Trying face " << ve.face << std::endl;
#endif
						if(faces_in_conflict_indices.find(ve.face) == faces_in_conflict_indices.end()){
							vertices[vj].face = ve.face;
							std::cout << "         Setting to face " << ve.face << std::endl;
							break;
						}
					}while(ve.face != old_face); // should never get back to first face
#ifdef PERIODIC_TRIANGULATION_DEBUG
					if(ve.face == old_face){
						std::cout << "      FAILED to find new face." << std::endl;
					}
#endif
				}
#ifdef PERIODIC_TRIANGULATION_DEBUG
				else{
					std::cout << "      Vertex " << vj << " ok with face: " << old_face << std::endl;
				}
#endif
			}
		}
		
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Hole: ";
		for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
			size_t v1 = faces[i->edge.face].v[(i->edge.across+2)%3];
			size_t v2 = faces[i->edge.face].v[(i->edge.across+1)%3];
			std::cout << "(" << v1 << "-" << v2 << ")o(" << i->offset[0] << "," << i->offset[1] << ") ";
		}std::cout << std::endl;
#endif
		// Remove all those faces.
		{
			// Update all indices
			// Find the mapping from old face indices to new indices
			std::vector<size_t> new_face_indices(faces.size());
			typename std::set<ParticularFace>::const_iterator cf = faces_in_conflict.begin();
			for(size_t i = 0, j = 0; i < faces.size(); ++i){
				if(cf != faces_in_conflict.end() && i > cf->face){
					++j;
					++cf;
				}
				new_face_indices[i] = i-j;
			}
#ifdef PERIODIC_TRIANGULATION_DEBUG
			for(size_t i = 0; i < faces.size(); ++i){
				std::cout << "   " << i << " -> " << new_face_indices[i] << std::endl;
			}
#endif
			// Remove all the faces
			{
				typename std::set<ParticularFace>::const_iterator i;
				size_t j = 0;
				for(i = faces_in_conflict.begin(); i != faces_in_conflict.end(); ++i, ++j){
					faces.erase(faces.begin()+(i->face)-j);
				}
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
				if(invalid_face != i->edge.face){
					i->edge.face = new_face_indices[i->edge.face];
				}
			}
			face_info.face = new_face_indices[face_info.face];
		}

#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Hole after face removal: ";
		for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
			size_t v1 = faces[i->edge.face].v[(i->edge.across+2)%3];
			size_t v2 = faces[i->edge.face].v[(i->edge.across+1)%3];
			std::cout << "(" << v1 << "-" << v2 << ")o(" << i->offset[0] << "," << i->offset[1] << ") ";
		}std::cout << std::endl;
#endif
		
		size_t v = vertices.size();
		vertices.push_back(nv);
		
		Offset nv_offset(-face_info.offset[0], -face_info.offset[1]);
		// nv_offset should NOT have negative components
		unsigned char nv_wrap = 0;
		if(nv_offset[0] > covering_extent[0]){ nv_wrap |= 1; }
		if(nv_offset[1] > covering_extent[1]){ nv_wrap |= 2; }
		// We can think of nv_wrap as whether or not the face that contained the vertex need wrapping,
		// or we can think of it as whether the vertex point need wrapping backwards.
		
		// Remarkably, connecting the new vertex to the boundary
		// vertices of the hole completes the triangulation.
		size_t first_new_face_index = faces.size();
		size_t num_new_faces = hole.size();
		size_t count = 0;
		faces.reserve(first_new_face_index + num_new_faces);
		for(typename Hole::iterator i = hole.begin(); i != hole.end(); ++i, ++count){
			size_t v1 = faces[i->edge.face].v[(i->edge.across+2)%3];
			unsigned char v1_wrap = faces[i->edge.face].GetVertexWrap((i->edge.across+2)%3);
			size_t v2 = faces[i->edge.face].v[(i->edge.across+1)%3];
			unsigned char v2_wrap = faces[i->edge.face].GetVertexWrap((i->edge.across+1)%3);
			unsigned char v_wrap = nv_wrap;
			unsigned char e_wrap = 0;
			if(i->offset[0]){ e_wrap |= 1; } // if offset is nonzero, it should only be negative (-2*covering_extent-1)
			if(i->offset[1]){ e_wrap |= 2; }
			// e_wrap is whether the neighboring face of the edge needed wrapping to the left/bottom
			// v1_wrap and v2_wrap are the wrapping bits for the vertices of the face of the current edge
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "      Adding with edge "
				<< v1 << "-" << v2 << ", wraps: " << int(nv_wrap) << int(e_wrap) << int(v1_wrap) << int(v2_wrap) << std::endl;
#endif
			faces.push_back(Face(
				v, v1, v2,
				i->edge.face, first_new_face_index+((count+1)%num_new_faces), first_new_face_index+((count+num_new_faces-1)%num_new_faces)
				));
			// Compute the offset information
			faces.back().wraps = MakeWrap(nv_wrap, e_wrap, v1_wrap, v2_wrap);
			
			// Patch up neighbors
			faces[i->edge.face].n[i->edge.across] = first_new_face_index+count;
		}
		vertices.back().face = first_new_face_index;

#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Done!" << std::endl;
#endif

		return v;
	}
	Edge EquivalentEdge(const Edge &e) const{ // given an edge, return an equivalent representation of the edge using the face opposite the edge
		size_t other_face = faces[e.face].n[e.across];
		int offset = faces[other_face].VertexOffset(faces[e.face].v[(e.across+1)%3]);
		return Edge(other_face, offset+1/*mod3 handled in constructor*/);
	}
	// Same as EquivalentEdge, but also provide it the offset of e.face, and this routine
	// will find the equivalent edge in e2, and also given the new offset of e2.face in e_face_offset
	// Note that e_face_offset can only change in units of 2*covering_extent at a time in each direction
	ParticularEdge EquivalentEdge(const ParticularEdge &e) const{
		size_t other_face = faces[e.edge.face].n[e.edge.across];
		int old_ev1 = (e.edge.across+1)%3;
		int old_ev2 = (e.edge.across+2)%3;
		int new_ev1 = faces[other_face].VertexOffset(faces[e.edge.face].v[old_ev1]);
		int new_ev2 = (new_ev1 + 2)%3;
		// new_ev1, new_ev2, new_across are the offsets of the vertices of the new face
		
		// If the wrapping status of the vertices of the edge are not the same in e and e2,
		// then e_face_offset will change.
		// If a wrapping bit changes from 1 to 0, then e_face_offset will increase by one unit
		// If a wrapping bit changes from 0 to 1, then e_face_offset will decrease by one unit
		
		
		// If one bit of the wrapping status of one vertex changes,
		// the other vertex's wrapping status bit should change in the same way
		/*
		unsigned char wraps[4];
		wraps[0] = faces[e.face].GetVertexWrap(old_ev1);
		wraps[1] = faces[e.face].GetVertexWrap(old_ev2);
		wraps[2] = faces[e2.face].GetVertexWrap(new_ev1);
		wraps[3] = faces[e2.face].GetVertexWrap(new_ev2);
		if(wraps[0] == wraps[2] && wraps[1] == wraps[3]){ return; }
		*/
		unsigned char wraps[2];
		wraps[0] = faces[e.edge.face].GetVertexWrap(old_ev1);
		wraps[1] = faces[other_face].GetVertexWrap(new_ev1);
		Offset e_face_offset(e.offset);
		if(wraps[0] != wraps[1]){
			for(int j = 1; j <= 2; ++j){
				if((0 == (wraps[0]&j)) && (0 != (wraps[1]&j))){
					e_face_offset[j-1] -= 2*covering_extent[j-1]+1;
				}else if((0 != (wraps[0]&j)) && (0 == (wraps[1]&j))){
					e_face_offset[j-1] += 2*covering_extent[j-1]+1;
				}
			}
		}
		return ParticularEdge(other_face, new_ev1+1/*mod3 handled in constructor*/, e_face_offset);
	}
	
	bool RemoveSingle(size_t vertex_index); // remove a single vertex (does not deal with copies)
	void MakeHole(size_t vertex_in_hole, Hole &hole){
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "In MakeHole, removing faces around vertex " << vertex_in_hole << ", adjacent face " << vertices[vertex_in_hole].face << std::endl;
#endif
		// Note that this function does not delete the vertex; it only removes all faces around a vertex
		std::vector<size_t> faces_in_hole;
		
		size_t cur_face = vertices[vertex_in_hole].face;
		size_t start = cur_face;
		ParticularEdge e(cur_face, (faces[cur_face].VertexOffset(vertex_in_hole)+2)%3);
		do{
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "      On face " << cur_face << std::endl;
#endif
			std::cout << "      Iteration edge between vertices " << faces[e.edge.face].v[(e.edge.across+1)%3] << "," << faces[e.edge.face].v[(e.edge.across+2)%3] << std::endl;
			// Move the hole boundary vertices' face away from the faces that will be deleted
			ParticularEdge en(EquivalentEdge(e.NextEdge()));
			//std::cout << "Neighbor edge between vertices " << faces[en.edge.face].v[(en.edge.across+1)%3] << "," << faces[en.edge.face].v[(en.edge.across+2)%3] << std::endl;
			size_t v;
			v = faces[cur_face].v[(en.edge.across+1)%3];
			if(vertices[v].face == cur_face){ vertices[v].face = en.edge.face; }
			v = faces[cur_face].v[(en.edge.across+2)%3];
			if(vertices[v].face == cur_face){ vertices[v].face = en.edge.face; }
			
			// Change the neighbors to indicate a deletion
			faces[en.edge.face].n[en.edge.across] = invalid_face;
			
			hole.push_back(en);
			faces_in_hole.push_back(cur_face);
			
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "   Hole:";
			for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
				std::cout << " " << (int)faces[i->edge.face].v[(i->edge.across+2)%3] << "-" << (int)faces[i->edge.face].v[(i->edge.across+1)%3] << "f" << i->edge.face;
			}std::cout << std::endl;
#endif
			e = EquivalentEdge(e.PrevEdge());
			cur_face = e.edge.face;
		}while(cur_face != start);
		
		std::sort(faces_in_hole.begin(), faces_in_hole.end());

#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Faces to remove:";
		for(size_t i = 0; i < faces_in_hole.size(); ++i){
			std::cout << " " << faces_in_hole[i];
		}std::cout << std::endl;
#endif

		// Update all indices
		// Find the mapping from old face indices to new indices
		std::vector<size_t> new_face_indices(faces.size());
		for(size_t i = 0, j = 0; i < faces.size(); ++i){
			if(j < faces_in_hole.size() && i > faces_in_hole[j]){
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
			if(invalid_face != i->edge.face){
				i->edge.face = new_face_indices[i->edge.face];
			}
		}
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   Hole:";
		for(typename Hole::const_iterator i = hole.begin(); i != hole.end(); ++i){
			std::cout << " " << (int)faces[i->edge.face].v[(i->edge.across+2)%3] << "-" << (int)faces[i->edge.face].v[(i->edge.across+1)%3] << "f" << i->edge.face;
		}std::cout << std::endl;
#endif
	}
	
	//// Specifically Delaunay primitives
	void FillHole(const Hole &hole_to_fill){
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "In FillHole, hole:";
		for(typename Hole::const_iterator i = hole_to_fill.begin(); i != hole_to_fill.end(); ++i){
			std::cout << " " << (int)faces[i->edge.face].v[(i->edge.across+2)%3] << "-" << (int)faces[i->edge.face].v[(i->edge.across+1)%3] << "o(" << i->offset[0] << "," << i->offset[1] << ")";
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
				for(size_t i = 0; i < 3; ++i){
					n[i] = h->edge.face;
					v[i] = faces[h->edge.face].v[(h->edge.across+2)%3];
					faces[h->edge.face].n[h->edge.across] = new_face_index; // update neighbors preemptively
					++h;
				}
				faces.push_back(Face(v[2],v[0],v[1], n[1],n[2],n[0]));
				// We consider this as joining v[2] with first edge of the hole
				unsigned char e_wrap = 0;
				if(hole.begin()->offset[0]){ e_wrap |= 1; }
				if(hole.begin()->offset[1]){ e_wrap |= 2; }
				unsigned char v1_wrap = faces[hole.begin()->edge.face].GetVertexWrap((hole.begin()->edge.across+2)%3);
				unsigned char v2_wrap = faces[hole.begin()->edge.face].GetVertexWrap((hole.begin()->edge.across+1)%3);
#ifdef PERIODIC_TRIANGULATION_DEBUG
				std::cout << "   Trivial filling: " << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
				std::cout << "      wraps: " << 0 << int(e_wrap) << int(v1_wrap) << int(v2_wrap) << std::endl;
#endif
				faces.back().wraps = MakeWrap(0, e_wrap, v1_wrap, v2_wrap);
				continue;
			}
			
			// We will now add a face with one of its edges being e
			ParticularEdge e = hole.front(); hole.pop_front();
			unsigned char e_wrap = 0;
			if(e.offset[0]){ e_wrap |= 1; }
			if(e.offset[1]){ e_wrap |= 2; }
			
			size_t v0 = faces[e.edge.face].v[(e.edge.across+2)%3];
			const unsigned char v0_wrap = faces[e.edge.face].GetVertexWrap((e.edge.across+2)%3);
			Pt2 p0 = GetFacePoint(e.GetParticularFace(), (e.edge.across+2)%3);
			size_t v1 = faces[e.edge.face].v[(e.edge.across+1)%3];
			const unsigned char v1_wrap = faces[e.edge.face].GetVertexWrap((e.edge.across+1)%3);
			Pt2 p1 = GetFacePoint(e.GetParticularFace(), (e.edge.across+1)%3);
			// The new triangle will be v0,v1,v2, where v2 is to be determined
			size_t v2;
			Pt2 p2;
			
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "   v0: " << v0 << ", v1: " << v1 << std::endl;
#endif
			
			// For all the other hole edges, we must now find the third vertex to join with it
			// We iterate over all the possible third vertices
			typename Hole::iterator hend = hole.end(); --hend; // the last vertex is on the edge one before the edge before e
			typename Hole::iterator h = hole.begin(); // hole.begin() is actually the edge after e since we popped away e
			typename Hole::iterator cut_location(h); // we will be cutting hole after this iterator location
			bool first = true;
			while(h != hend){
				int which = (h->edge.across+1)%3;
				size_t cur_vertex = faces[h->edge.face].v[which];
				// cur_vertex will take on all the proper possible values for v2 along the hole boundary
				// We will go through all of them without early exit to get the last one
				
				// Throughout this iteration, v2 will store the last (in
				// the sense of previous) possible vertex to use.
				
				// Now we only have to check the Delaunay criterion.
				Pt2 p = GetFacePoint(h->GetParticularFace(), which);
#ifdef PERIODIC_TRIANGULATION_DEBUG
				std::cout << "      Testing Delaunay for vertex " << cur_vertex << " at (" << p.r[0] << "," << p.r[1] << ")" << std::endl;
#endif
				if(Orient2(p0, p1, p) > 0){ // if the proposed face has the proper orientation...
					if(first || InCircle2(p0, p1, p2, p) > 0){ // If p is INSIDE the last circumcircle, then we better use p instead!
						v2 = cur_vertex;
						p2 = p;
						cut_location = h;
						first = false;
					}
				}
				
				++h;
			}
	
			unsigned char v2_wrap = 0;
			{
				Offset v2_offset(e.GetParticularFace().offset);
				// v2_offset should NOT have negative components
				if(v2_offset[0] > covering_extent[0]){ v2_wrap |= 1; }
				if(v2_offset[1] > covering_extent[1]){ v2_wrap |= 2; }
			}
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "   Adding face with vertices " << v0 << "," << v1 << "," << v2 << " (not in this order) Wraps: " << int(v2_wrap) << int(e_wrap) << int(v0_wrap) << int(v1_wrap) << std::endl;
#endif
			
			// We now need to make a face defined by edge e and vertex v2
			// We also need to split the hole into two holes since the triangle could have
			// not clipped an ear off the hole.
			
			// First we check for the two cases where the hole is not actually split.
			ParticularEdge en(e); // the neighboring edge directly after/before e (we have to initialize it with something, so use e)
			int offset;
			if( // assign en to edge directly after e (in CCW order) and see if the far vertex is v2
				((en = hole.front()), (v2 == faces[en.edge.face].v[(en.edge.across+1)%3]))
			){ // If v2 is directly after v0 and v1, then there is no hole split
#ifdef PERIODIC_TRIANGULATION_DEBUG
				std::cout << "   In case 1; no split, new face: " << v2 << ", " << v0 << ", " << v1 << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v2,v0,v1, e.edge.face,en.edge.face,invalid_face));
				faces.back().wraps = MakeWrap(v2_wrap, e_wrap, v0_wrap, v1_wrap);

				// Update neighbors
				faces[e.edge.face].n[e.edge.across] = new_face_index;
				faces[en.edge.face].n[en.edge.across] = new_face_index;
				
				Offset new_offset = e.offset;
				for(size_t j = 0; j < 2; ++j){
					if(new_offset[j] < en.offset[j]){ new_offset[j] = en.offset[j]; }
				}
				
				// The hole's old edge needs replacing with the new face's edge
				hole.pop_front();
				hole.push_front(ParticularEdge(new_face_index, 2, new_offset));
				holes.push_front(hole);
			}else if( // assign en to edge directly before e (in CCW order) and see if the far vertex is v2
				((en = hole.back()), (v2 == faces[en.edge.face].v[(en.edge.across+2)%3]))
			){ // If v2 is directly before v0 and v1, then there is no hole split
#ifdef PERIODIC_TRIANGULATION_DEBUG
				std::cout << "   In case 2; no split, new face: " << v2 << ", " << v0 << ", " << v1 << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v2,v0,v1, e.edge.face, invalid_face,en.edge.face));
				faces.back().wraps = MakeWrap(v2_wrap, e_wrap, v0_wrap, v1_wrap);
				// Update neighbors
				faces[e.edge.face].n[e.edge.across] = new_face_index;
				faces[en.edge.face].n[en.edge.across] = new_face_index;
				
				Offset new_offset = e.offset;
				for(size_t j = 0; j < 2; ++j){
					if(new_offset[j] < en.offset[j]){ new_offset[j] = en.offset[j]; }
				}
				
				// The hole's old edge needs replacing with the new face's edge
				hole.pop_back();
				hole.push_front(ParticularEdge(new_face_index, 1, new_offset));
				holes.push_front(hole);
			}else{ // We will have a hole split
#ifdef PERIODIC_TRIANGULATION_DEBUG
				std::cout << "   Splitting hole" << std::endl;
#endif
				size_t new_face_index = faces.size();
				// Create the new face
				faces.push_back(Face(v0,v1,v2, invalid_face,invalid_face,e.edge.face));
				faces[e.edge.face].n[e.edge.across] = new_face_index; // update neighbor
				
				Offset new_offset = e.offset;
				
				// Split the hole
				Hole new_hole;
				for(size_t j = 0; j < 2; ++j){
					if(new_offset[j] < cut_location->offset[j]){ new_offset[j] = cut_location->offset[j]; }
				}
				++cut_location;
				for(size_t j = 0; j < 2; ++j){
					if(new_offset[j] < cut_location->offset[j]){ new_offset[j] = cut_location->offset[j]; }
				}
				// Not sure if the computation of new_offset here is correct
				
				{
					v2_wrap = 0;
					if(new_offset[0]){ v2_wrap |= 1; }
					if(new_offset[1]){ v2_wrap |= 2; }
				}
				unsigned char wraps = MakeWrap(v2_wrap, e_wrap, v0_wrap, v1_wrap);
				faces.back().wraps = (((wraps & 0x3C) >> 2) | ((wraps & 0x03) << 4));
				
				while(hole.begin() != cut_location){
					new_hole.push_back(hole.front()); // note this has to be push_back to keep the orientation the same
					hole.pop_front();
				}
				hole.push_front(ParticularEdge(new_face_index, 1, new_offset));
				new_hole.push_front(ParticularEdge(new_face_index, 0, new_offset));
				holes.push_front(hole);
				holes.push_front(new_hole);
			}
		}
	}
	

	NumericType InCircumcircle(const ParticularFace &face_info, const LatticePoint &p){
		Pt2 fp[3];
		GetFacePoints(face_info, fp);
		return InCircle2(
			fp[0], fp[1], fp[2],
			L(
				p.st[0]+NumericType(p.offset[0]),
				p.st[1]+NumericType(p.offset[1]))
			);
	}
	NumericType Orientation(size_t v0, size_t v1, const LatticePoint &p){
		return Orient2(
			L(
				vertices[v0].p.st[0]+NumericType(vertices[v0].p.offset[0]),
				vertices[v0].p.st[1]+NumericType(vertices[v0].p.offset[1])),
			L(
				vertices[v1].p.st[0]+NumericType(vertices[v1].p.offset[0]),
				vertices[v1].p.st[1]+NumericType(vertices[v1].p.offset[1])),
			L(
				p.st[0]+NumericType(p.offset[0]),
				p.st[1]+NumericType(p.offset[1]))
			);
	}
	unsigned char MakeWrap(unsigned char nv_wrap, unsigned char e_wrap, unsigned char v1_wrap, unsigned char v2_wrap) const{
		// This function computes the wraps parameter for a newly created face.
		// Each parameter is a 2 bit field, the low order bit being the wrap in the [0] index direction
		// and the second bit being the wrap in the [1] direction.
		//
		// It is assumed the face was created by joining a vertex with an edge.
		// The edge belongs to a face on the (OTHER side of the edge) than the vertex we are joining with.
		// The edge must be a particular edge (an edge with an offset).
		// nv_wrap is the vertex's wrap relative to the canonical face positions (faces with zero offset).
		// Therefore nv_wrap should be cleared if the vertex lies in the fundamental parallelogram within a face,
		// and set if the vertex lies on the far left or bottom within the fundamental parallelogram,
		// necessetating faces at the top/right to be shifted backwards to contain the vertex.
		//
		// e_wrap is set if the particular edge's offset is nonzero. Note that the edge's offset
		// can only be zero or -2*covering_extent-1. It cannot be greater than zero or less than
		// -2*convering_extent-1 or else the face would lie entirely outside of the fundamental parallelogram.
		// Furthermore, it must be a multiple of 2*covering_extent+1 obviously since it wraps around
		// the covering domain.
		//
		// v1_wrap and v2_wrap are the (vertex offsets pulled from the face of the particular edge)
		// of the two vertices at the ends of the edge. They should be obtained by
		// faces[cur_partic_edge.edge.face].GetVertexWrap(i)
	
		// Truth table for wrapping:
		// nv_wrap
		// | e_wrap
		// | | v1_wrap
		// | | | v2_wrap
		// | | | | desired final wraps bitfield (in order of v, v1, v2)
		// | | | | |
		// 1 1 0 0 100 A1
		// 1 1 0 1 101 A2
		// 1 1 1 0 110 A3
		// 1 0 0 0 000 A4
		// 0 1 0 0 100 B1
		// 0 1 0 1 101 B2
		// 0 1 1 0 110 B3
		// 0 0 0 0 000 B4
		// 0 0 0 0 000 C1
		// 0 0 0 1 001 C2
		// 0 0 1 0 010 C3
		// 0 1 0 0 011 C4
		// 0 1 1 1 000 D
		//
		// The table is simple; the final wrapping bits are always just
		// the concatenation of e_wrap, v1_wrap, and v2_wrap, EXCEPT
		// for line C4 and D
		
		// These cases are visualized as follows:
		//
		//         +-------------------------+
		//         |                         |
		//        2|                        2|
		//     ,--===--.                 ,--===--.
		//   ,'    |    `.             ,'    |    `.
		// ||    A | B    ||         ||    C |      ||
		// ||1   + | +    ||4       1||    + |      ||4
		// ||      |      ||         ||      |      ||
		//   `.    |    ,'             `.    |    ,'
		//     '--===--'                 '--===--'
		//        3|                        3|
		//         |                         |
		//         |                         |
		//         +-------------------------+
		// We only look at one direction here; the two directions can be treated separately.
		// The letter indicates the position of the vertex nv (determined by nv_wrap)
		// The numbered positions indicate the position of an edge of the hole we will
		// connect in the loop below, as well as the wrap parameters of v1 and v2.
		// Case D is special in that we should never have a face where all vertices are wrapped
		
		unsigned char wraps = (e_wrap | (v1_wrap << 2) | (v2_wrap << 4));
		if(
			((wraps & 0x15) == 0x15) // special case D
			|| (!(nv_wrap&1) && ((wraps & 0x01) == 0x01)) // special case C4
			){ wraps ^= 0x15; }
		if(
			((wraps & 0x2A) == 0x2A)
			|| (!(nv_wrap&2) && ((wraps & 0x02) == 0x02))
			){ wraps ^= 0x2A; }
		return wraps;
	}
private:
	// Embedded GameRand implementation for location queries
	mutable unsigned int rng_low, rng_high;
	unsigned int rng() const{
		rng_high = (rng_high << 16) + (rng_high >> 16);
		rng_high += rng_low; rng_low += rng_high;
		return rng_high;
	}
};


// At all times we must have at least 3 vertices in the triangulation
template <typename NumericType>
Triangulation2<NumericType>::Triangulation2(const Lattice &lattice, std::vector<Vertex> &existing_vertices):
	L(lattice),
	one_sheeted(false)
{
	// initialize GameRand
	rng_high = 0xDEADBEEF;
	rng_low = rng_high ^ 0x49616E42;
	
	if(existing_vertices.size() < 1){ return; }


	// Find the unit cell offsets need to guarantee a simplicial complex triangulation
	// For a generic lattice, the covering set is not just the 3 unit cells in each direction
	// We need to go as far out as 3 times the longest diagonal length of the fundamental parallelogram
	covering_extent[0] = 0; covering_extent[1] = 0;
	{ // scoping
		// In order to robustly find end of the covering set, we don't just
		// measure distance from the origin. We know that the longest diagonal
		// of the fundamental parallelogram will be u+v or u-v.
		Coord2 final_coord;
		if((L(-1,1)-Pt2::Origin).LengthSq() < (L(1,1)-Pt2::Origin).LengthSq()){
			final_coord = Coord2(1,1);
		}else{
			final_coord = Coord2(-1,1);
		}
		
		int st[2] = {0,0};
		NumericType r2 = 0;
		bool final_found = false;
		do{
			std::vector<Coord2> ring_pts;
			L.GetNextPointRing(st[0], st[1], ring_pts);
			for(typename std::vector<Coord2>::const_iterator i = ring_pts.begin(); i != ring_pts.end(); ++i){
				for(int j = 0; j < 2; ++j){
					int a = abs(i->v[j]);
					if(covering_extent[j] < a){ covering_extent[j] = a; }
				}
				if(*i == final_coord){ final_found = true; }
			}
		}while(!final_found);
	}
	// The number of sheets in each direction for the covering space is
	// 2*covering_extent[i]+1, and we will use 0-centric indexing, so
	// almost half the offsets will be negative numbers.
	
	// Set up the initial triangulation
	// We will need to link one point in the fundamental parallelogram to
	// one point in each of 3 parallelograms in the 1-ring around the
	// fundamental one. We will refer to these parallelograms by their
	// offsets; the fundamental being (0,0).
	// Because the lattice basis vectors are reduced, we always have
	// the links (0,0)->(1,0) and (0,0)->(0,1). The third is determined
	// by the lengths of the diagonals across the fundamental
	// parallelogram.
	Vec2 u,v;
	L.GetBasis(u,v);
	NumericType upv2 = (u+v).LengthSq(), umv2 = (u-v).LengthSq();
	bool link_umv = false;
	if(umv2 < upv2){
		link_umv = true;
	}
	
	

	// The edge length threshold is such that if all edges are shorter than the threshold,
	// it is guaranteed that there will be no self edges.
	// Following the Caroli & Teillaud paper, assume there is an empty circle C of diameter d.
	// The longest edge of a face is bounded below by the edge length of a regular triangle
	// with circumscribing circle C, which is d*sqrt(3)/2. We need the diameter of any empty
	// circle to be smaller than D/2, where D is the diameter of the largest inscribed circle
	// within the fundamental parallelogram.
	// The threshold^2 is then D^2*3/16
	
	// Let a be the length of the shorter of the two basis vectors, let p = 0.5|u+v|, m = 0.5|u-v|,
	// and r is the radius of D that we are after:
	//             _
	//         +  |\          ---
	//        /|\   \         /|\
	//       / | \   \l        |
	//      /  |  \   \        |
	//    a/   |   \   \       | p
	//    / ___+___ \  -\|     |
	//   /,'   |   ',\         |
	//  //     |     \\       \|/
	// +-|-----+-----|-+      ---
	//   |<-r->|<--m-->|
	//
	// Let l be the distance from the corner of the parallelogram (rhombus, it's all the same)
	// to the point of tangency with the inscribed circle.
	// Then we have
	//   l^2 + r^2 = p^2        = a^2 - m^2
	//   (a-l)^2 + r^2 = m^2
	// Eliminating l^2,
	//   l^2 + r^2 = a^2 - m^2
	//   a^2 + l^2 - 2al + r^2 = m^2
	//   -->  a^2 - al = m^2
	//        l = (a^2-m^2)/a = p^2/a
	// From the very first equation,
	//   r^2 = p^2 - (p^2)^2/a^2 = (a^2 - m^2) - (a^4 - 2a^2m^2 + m^4)/a^2 = m^2 (1 - m^2/a^2)
	{
		// This is coded in a way that is closed under arithmetic on rational_radical's
		NumericType m2 = (link_umv ? umv2 : upv2); m2 /= NumericType(4);
		NumericType a2 = u.LengthSq(), v2 = v.LengthSq();
		if(v2 < a2){ a2 = v2; }
		NumericType D2_4 = m2*(NumericType(1) - m2/a2); // The real D^2 is 4 times of this
		edge_length_threshold2 = D2_4 * (NumericType(3)/NumericType(4));
	}
	
#ifdef PERIODIC_TRIANGULATION_DEBUG
	std::cout << "Covering extent: " << covering_extent[0] << "," << covering_extent[1] << std::endl;
#endif

	Offset st;
	int nst[2];
	for(int i = 0; i < 2; ++i){ nst[i] = 2*covering_extent[i]+1; }
	int num_copies = nst[0]*nst[1];
	
	vertices.reserve(size_t(num_copies) * existing_vertices.size());
	// Insert all copies of first vertex
	size_t center_vertex_index = covering_extent[1]*(2*covering_extent[0]+1)+covering_extent[0];
	vertex_to_copies[center_vertex_index].reserve(num_copies);
	for(st[1] = -covering_extent[0]; st[1] <= covering_extent[0]; ++st[1]){
		for(st[0] = -covering_extent[0]; st[0] <= covering_extent[0]; ++st[0]){
			size_t new_vertex_index = vertices.size();
			size_t face_index = 2*new_vertex_index;
			copy_to_vertex[new_vertex_index] = center_vertex_index;
			vertex_to_copies[center_vertex_index].push_back(new_vertex_index);
			vertices.push_back(Vertex(existing_vertices.front(), st));
			vertices.back().face = face_index;
		}
	}
	for(st[1] = -covering_extent[0]; st[1] <= covering_extent[0]; ++st[1]){
		for(st[0] = -covering_extent[0]; st[0] <= covering_extent[0]; ++st[0]){
			size_t v0  =  (covering_extent[1]+st[1]+0)*nst[0]+(covering_extent[0]+st[0]+0);
			size_t v1  =  (covering_extent[1]+st[1]+0)*nst[0]+(covering_extent[0]+st[0]+1)%nst[0];
			size_t v1n =  (covering_extent[1]+st[1]+0)*nst[0]+(covering_extent[0]+st[0]+nst[0]-1)%nst[0];
			size_t v2  = ((covering_extent[1]+st[1]+1)%nst[1])*nst[0]+(covering_extent[0]+st[0]+0);
			size_t v2n = ((covering_extent[1]+st[1]+nst[1]-1)%nst[1])*nst[0]+(covering_extent[0]+st[0]+0);
			size_t v3  = ((covering_extent[1]+st[1]+1)%nst[1])*nst[0]+(covering_extent[0]+st[0]+1)%nst[0];
			size_t face_base = 2*faces.size();
			//unsigned char wrap0 = 0; // v0 is never wrapped
			unsigned char wrap1 = ((st[0] == covering_extent[0]) ? 0x1 : 0x0);
			unsigned char wrap2 = ((st[1] == covering_extent[1]) ? 0x2 : 0x0);
			unsigned char wrap3 = ((st[0] == covering_extent[0]) ? 0x1 : 0x0);
			              wrap3|= ((st[1] == covering_extent[1]) ? 0x2 : 0x0);
			if(link_umv){
				faces.push_back(Face(
					v0,v1,v2,
					2*v0+1, 2*v1n+1, 2*v2n+1,
					(wrap1 << 2) | (wrap2 << 4)
					));
				faces.push_back(Face(
					v1,v3,v2,
					2*v2, 2*v0, 2*v1+1,
					wrap1 | (wrap3 << 2) | (wrap2 << 4)
					));
			}else{
				faces.push_back(Face(
					v0,v1,v3,
					2*v1+1, 2*v0+1, 2*v2n+1,
					(wrap1 << 2) | (wrap3 << 4)
					));
				faces.push_back(Face(
					v0,v3,v2,
					2*v2, 2*v1n, 2*v0,
					(wrap3 << 2) | (wrap2 << 4)
					));
			}
		}
	}
	
	// Insert remaining faces
	AddVertices(existing_vertices.begin()+1,existing_vertices.end());
}

template <typename NumericType>
Triangulation2<NumericType>::Triangulation2(const Triangulation2& t):
	L(t.L),
	vertices(t.vertices),
	faces(t.faces),
	one_sheeted(t.one_sheeted),
	rng_low(t.rng_low),
	rng_high(t.rng_high)
{
}

template <typename NumericType>
typename Triangulation2<NumericType>::Pt2 Triangulation2<NumericType>::GetPoint(size_t vertex_index) const{
	return L(
		vertices[vertex_index].p.st[0]+NumericType(vertices[vertex_index].p.offset[0]),
		vertices[vertex_index].p.st[1]+NumericType(vertices[vertex_index].p.offset[1]));
}

template <typename NumericType>
void Triangulation2<NumericType>::GetFacePoints(const ParticularFace &face_info, Pt2 v[3]) const{
	// Faces can never wrap to the left/bottom

	for(size_t i = 0; i < 3; ++i){
		v[i] = GetFacePoint(face_info, i);
	}
}

template <typename NumericType>
typename Triangulation2<NumericType>::Pt2 Triangulation2<NumericType>::GetFacePoint(const ParticularFace &face_info, int which) const{
	size_t cur_vertex_index = faces[face_info.face].v[which];
	Offset cur_offset(vertices[cur_vertex_index].p.offset);
	for(size_t j = 1; j <= 2; ++j){
		if(0 != (faces[face_info.face].GetVertexWrap(which)&j)){
			cur_offset[j-1] += 2*covering_extent[j-1]+1;
		}
	}
	cur_offset += face_info.offset;
	return L(
		vertices[cur_vertex_index].p.st[0]+NumericType(cur_offset[0]),
		vertices[cur_vertex_index].p.st[1]+NumericType(cur_offset[1]));
}

template <typename NumericType>
typename Triangulation2<NumericType>::ParticularFace Triangulation2<NumericType>::ContainingFace(const LatticePoint &p, size_t face_guess) const{
	// Make sure the point is in the domain
	if(one_sheeted){
		if(!(p.offset == Offset()) || (p.st[0] < NumericType(0)) || (p.st[0] > NumericType(1))){ return ParticularFace(invalid_face); }
	}else{
		if(p.offset.off[0] < -covering_extent[0] || p.offset.off[0] > covering_extent[0]){ return ParticularFace(invalid_face); }
		if(p.offset.off[1] < -covering_extent[1] || p.offset.off[1] > covering_extent[1]){ return ParticularFace(invalid_face); }
		if((p.st[0] < NumericType(0)) || (p.st[0] > NumericType(1))){ return ParticularFace(invalid_face); }
	}
#ifdef PERIODIC_TRIANGULATION_DEBUG
	std::cout << "In ContainingFace, looking for " << p.offset[0] << "+" << p.st[0] << "," << p.offset[1] << "+" << p.st[1] << std::endl;
	int n_iter = 0;
#endif
	ParticularFace face_info(face_guess);
	while(true){
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "   In face " << face_info.face << " with offset " << face_info.offset[0] << "," << face_info.offset[1] << std::endl;
		++n_iter;
		if(n_iter > 10*faces.size()){ break; }
#endif
		int outside[2]; int noutside = 0;
		Pt2 fv[3];
		GetFacePoints(face_info, fv);
		for(size_t i = 0; i < 3; ++i){
			if(Orient2(fv[i], fv[(i+1)%3], L(p.st[0]+NumericType(p.offset[0]), p.st[1]+NumericType(p.offset[1]))) < 0){
				outside[noutside++] = (i+2)%3; //faces[containing_face].n[(i+2)%3];
				if(noutside >= 2){ break; }
			}
		}

#ifdef PERIODIC_TRIANGULATION_DEBUG
		for(size_t i = 0; i < noutside; ++i){
			std::cout << "      Point is outside edge across from vertex " << faces[face_info.face].v[outside[i]] << std::endl;
		}
#endif

		if(0 == noutside){
#ifdef PERIODIC_TRIANGULATION_DEBUG
			std::cout << "   Found!" << std::endl;
#endif
			break;
		}else{
			int iflip = 0;
			if(noutside > 1){
				// pick a random side to flip over to
				iflip = rng()&1;
			}
			//containing_face = outside[iflip];
			
			//replace this with EquivalentEdgeWrap
			ParticularEdge e2(EquivalentEdge(ParticularEdge(face_info, outside[iflip])));
			face_info.face = e2.edge.face;
			face_info.offset = e2.offset;
		}
	}
	return face_info;
}

template <typename NumericType>
size_t Triangulation2<NumericType>::AddVertex(const Vertex &v, size_t *face_index){
	size_t ret;
	size_t f = ((NULL == face_index) ? 0 : *face_index);
	if(one_sheeted){
		ParticularFace face_info(ContainingFace(v.p, f));
		if(NULL != face_index){ *face_index = f; } // set it back for feedback
		
		ret = AddVertexInFace(v, face_info);
	}else{
		Offset st;
		int nst[2];
		for(int i = 0; i < 2; ++i){ nst[i] = 2*covering_extent[i]+1; }
		
		// Insert all copies of first vertex
		bool first = true;
		size_t center_vertex_index = vertices.size() + covering_extent[1]*(2*covering_extent[0]+1) + covering_extent[0];
		vertex_to_copies[center_vertex_index].reserve(nst[0]*nst[1]);
		for(st[1] = -covering_extent[1]; st[1] <= covering_extent[1]; ++st[1]){
			for(st[0] = -covering_extent[0]; st[0] <= covering_extent[0]; ++st[0]){
				Offset f_offset;
				Vertex nv(v, st);
				ParticularFace face_info(ContainingFace(nv.p, f));
				if(NULL != face_index){ *face_index = face_info.face; } // set it back for feedback
				
				size_t idx = AddVertexInFace(Vertex(nv, st), face_info);
				vertex_to_copies[center_vertex_index].push_back(idx);
				copy_to_vertex[idx] = center_vertex_index;
				
				f = vertices[idx].face;
				if(first){ ret = idx; first = false; }
				/*
				if(st[1] == -covering_extent[1]+0 && st[0] == -covering_extent[0]+2){
					return ret;
				}
				//*/
			}
		}
	}
	return ret;
}

template <typename NumericType>
template <class InputIterator>
size_t Triangulation2<NumericType>::AddVertices(InputIterator begin, InputIterator end){
	size_t face_index = 0;
	size_t ret; bool ret_found = false;
	for(InputIterator i = begin; i != end; ++i){
		size_t idx = AddVertex(*i, &face_index);
		if(false == ret_found){ ret = idx; ret_found = true; }
	}
	return ret;
}

template <typename NumericType>
void Triangulation2<NumericType>::Flip(const Edge &e){
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
size_t Triangulation2<NumericType>::Split(const Edge &e, int tag){ // returns new vertex index
	Edge e2(EquivalentEdge(e));
	Pt2 p0(vertices[faces[e.face].v[(e.across+1)%3]].p);
	Pt2 p1(vertices[faces[e.face].v[(e.across+1)%3]].p);
	p0 += (NumericType(1)/NumericType(2)) * (p1-p0);
	size_t v = AddVertexInFace(Vertex(p0, tag), e.face);
	Flip(e2);
	return v;
}

template <typename NumericType>
bool Triangulation2<NumericType>::Remove(size_t vertex_index){
	if(vertices.size() < 1){ return false; }
	
	if(one_sheeted){
		return RemoveSingle(vertex_index);
	}else{
		// Find all copies of the vertex
		size_t center_vertex_index = copy_to_vertex[vertex_index];
		std::vector<size_t> vertices_to_remove(vertex_to_copies[center_vertex_index]);
		// Since vertex indices after the removed vertex will change after removal,
		// remove from largest to smallest index.
		std::sort(vertices_to_remove.begin(), vertices_to_remove.end());
#ifdef PERIODIC_TRIANGULATION_DEBUG
		std::cout << "Removing vertices:";
		for(size_t i = 0; i < vertices_to_remove.size(); ++i){
			std::cout << " " << vertices_to_remove[i];
		}
		std::cout << std::endl;
#endif
		RemoveSingle(34); return true;
		// const_reverse_iterator cannot work here (http://gcc.gnu.org/ml/gcc-bugs/1998-05/msg00738.html)
		for(std::vector<size_t>::reverse_iterator i = vertices_to_remove.rbegin(); i != vertices_to_remove.rend(); ++i){
			copy_to_vertex.erase(*i);
			RemoveSingle(*i);
			/*
			if(*i == 6){
				break;
			}
			//*/
		}
		vertex_to_copies.erase(center_vertex_index);
		
		return true;
	}
}

template <typename NumericType>
bool Triangulation2<NumericType>::RemoveSingle(size_t vertex_index){
	if(vertices.size() < 1){ return false; }
	Hole hole;
	MakeHole(vertex_index, hole);
	FillHole(hole);
	
	for(size_t i = 0; i < faces.size(); ++i){
		for(size_t j = 0; j < 3; ++j){
			if(faces[i].v[j] > vertex_index){ --faces[i].v[j]; }
		}
	}
	vertices.erase(vertices.begin()+vertex_index);
	return true;
}

}; // namespace Periodic

#endif // _PERIODIC_TRIANGULATION_HPP_
