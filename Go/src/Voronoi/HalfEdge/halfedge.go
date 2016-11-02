package halfedge

import (
    //"fmt"
)

type Vector struct {
    X   float32
    Y   float32
    //Z   float32
}

var InfiniteFace = &HEFace{}
var EmptyEdge    = &HEEdge{}
var EmptyVertex  = &HEVertex{}

// Implement the sort.Interface, so it can be sorted by the standard algorithm of Go.
type PointList []Vector
func (slice PointList) Len() int {
    return len(slice)
}
func (slice PointList) Less(i, j int) bool {
    return slice[i].X < slice[j].X;
}
func (slice PointList) Swap(i, j int) {
    slice[i], slice[j] = slice[j], slice[i]
}


type HEVertex struct {
    Pos             *Vector
    // It is possible, several edges leave from this vertex.
    // Then it points to any one of them, arbitrarily.
    ELeaving        *HEEdge
}

type HEFace struct {
    // Especially for Voronoi
    ReferencePoint  *Vector

    // Points to an arbitrary edge of its polygon
    EEdge           *HEEdge
}

type HEEdge struct {
    VOrigin         *HEVertex
    ETwin           *HEEdge
    // Starts from: this->ETwin->VOrigin. With this->FFace == this->ENext->FFace
    // Following ENext will traverse the polygon around FFace!
    ENext           *HEEdge
    FFace           *HEFace
}



// Todo but not important right now.
func (f *HEFace) CalcNormal() Vector {
    return Vector{}
}







