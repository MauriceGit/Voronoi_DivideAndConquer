package halfedge

import (
    //"fmt"
    . "Voronoi/Vector"
)

// Just to make the code more readable.
type EdgeIndex    int
type FaceIndex    int
type VertexIndex  int

var EmptyFace    = FaceIndex(-1)
var EmptyEdge    = EdgeIndex(-1)
var EmptyVertex  = VertexIndex(-1)


/**
 * All edges, vertices and faces are kept in static, not-resizable slices
 * and therefore tight data representation.
 *
 * References to vertices, edges or faces are kept as indices into the corresponding slices.
 *
 * It is guaranteed, that the initial lists are long enough for the to-be constructed voronoi.
 * Data does not need to be relocated for normal voronoi creation.
 */

type HEVertex struct {
    Pos             Vector
    // It is possible, several edges leave from this vertex.
    // Then it points to any one of them, arbitrarily.
    // This shit doesn't help anyone, does it?
    ELeaving        EdgeIndex
}

type HEFace struct {
    // Especially for Voronoi
    ReferencePoint  Vector

    // Points to an arbitrary edge of its polygon
    // Only arbitrary for closed faces!!!!
    // Otherwise it HAS to point to the FIRST edge (counter clockwise)!!!
    EEdge           EdgeIndex

    // For convenience, so we don't have to recalculate all the time.
    //IsClosed        bool
}

type HEEdge struct {

    VOrigin         VertexIndex
    ETwin           EdgeIndex
    // Starts from: this->ETwin->VOrigin. With this->FFace == this->ENext->FFace
    // Following ENext will traverse the polygon around FFace!
    ENext           EdgeIndex
    FFace           FaceIndex

    // Iff VOrigin is NOT defined, an edge representation is required
    // (general direction + some kind of starting point that is NOT an official HEVertex!)
    TmpEdge         Edge
}

// Todo but not important right now.
func (f *HEFace) CalcNormal() Vector {
    return Vector{}
}








