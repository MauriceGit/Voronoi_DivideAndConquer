package halfedge

import (
    //"fmt"
    . "Voronoi/Vector"
)

var EmptyFace    = &HEFace{}
var EmptyEdge    = &HEEdge{}
var EmptyVertex  = &HEVertex{}

type HEVertex struct {
    Pos             Vector
    // It is possible, several edges leave from this vertex.
    // Then it points to any one of them, arbitrarily.
    // This shit doesn't help anyone, does it?
    ELeaving        *HEEdge
}

type HEFace struct {
    // Especially for Voronoi
    ReferencePoint  Vector

    // Points to an arbitrary edge of its polygon
    // Only arbitrary for closed faces!!!!
    // Otherwise it HAS to point to the FIRST edge (counter clockwise)!!!
    EEdge           *HEEdge

    // For convenience, so we don't have to recalculate all the time.
    IsClosed        bool
}

type HEEdge struct {

    VOrigin         *HEVertex
    ETwin           *HEEdge
    // Starts from: this->ETwin->VOrigin. With this->FFace == this->ENext->FFace
    // Following ENext will traverse the polygon around FFace!
    ENext           *HEEdge
    FFace           *HEFace

    // Iff VOrigin is NOT defined, an edge representation is required
    // (general direction + some kind of starting point that is NOT an official HEVertex!)
    TmpEdge         Edge
}

// Todo but not important right now.
func (f *HEFace) CalcNormal() Vector {
    return Vector{}
}

// Creates a line from the HEEdge. Depeding on its state from the edge or the TmpEdge.
func (e HEEdge) Line (amplified bool) Edge {
    if e.VOrigin == EmptyVertex || e.ETwin.VOrigin == EmptyVertex {
        if amplified {
            tmpE := e.TmpEdge.Copy()
            tmpE.Amplify(100.0)
            return tmpE
        } else {
            return e.TmpEdge
        }
    } else {
        return Edge {
            Pos: e.VOrigin.Pos,
            Dir: Sub(e.ETwin.VOrigin.Pos, e.VOrigin.Pos),
        }
    }
}






