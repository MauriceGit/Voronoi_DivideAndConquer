package main

import (
    . "Voronoi/Vector"
    . "Voronoi/HalfEdge"
    "fmt"
    "math/rand"
    "sort"
    //"math"
)

type Voronoi struct {
    heVertexList    []HEVertex
    heEdgeList      []HEEdge
    heFaceList      []HEFace
}

type ConvexHull []*HEFace

func (v Voronoi) pprintShort () {
    fmt.Printf("Voronoi: ")
    fmt.Printf("   Vertices: %v\n", len(v.heVertexList))
    fmt.Printf("   Edges   : %v\n", len(v.heEdgeList))
    fmt.Printf("   Faces   : %v\n", len(v.heFaceList))
}
func (v Voronoi) pprint () {
    fmt.Printf("Voronoi: ")
    fmt.Printf("   Vertices: %v\n", v.heVertexList)
    fmt.Printf("   Edges   : %v\n", v.heEdgeList)
    fmt.Printf("   Faces   : %v\n", v.heFaceList)
}

// Calculates a convex hull for the given voronoi. Runs in O(n).
// The convex hull consists of a list of Reference points.
// As there is an exactly 1:1 representation of Reference points and HEFaces,
// a list of HEFaces is returned in Order (clockwise, as edges are counter-clockwise) of the convex hull.
func (v Voronoi)ConvexHull() ConvexHull {
    var outerFace = InfiniteFace

    // find outer face as starting point for convex hull.
    for i,face := range v.heFaceList {
        if !face.IsClosed {
            outerFace = &v.heFaceList[i]
            break
        }
    }
    var firstFace = outerFace
    var convexHullList = []*HEFace{outerFace}

    lastEdge := outerFace.EEdge

    if lastEdge != EmptyEdge {

        // find last edge of that face.
        for *lastEdge.ENext != *EmptyEdge {
            lastEdge = lastEdge.ENext
        }

        nextFace := lastEdge.ETwin.FFace

        // Now lastEdge is either EmptyEdge or the last edge before infinity regarding that face.
        for *nextFace != *firstFace {
            convexHullList = append(convexHullList, nextFace)
            lastEdge = nextFace.EEdge
            // find last edge of that face.
            for lastEdge.ENext != EmptyEdge {
                lastEdge = lastEdge.ENext
            }
            nextFace = lastEdge.ETwin.FFace
        }
    }

    return convexHullList
}

// There can only be 1 or 2 points in the list!!! This HAS to be guaranteed!
func createTrivialVoronoi(points PointList) Voronoi {

    if len(points) == 1 {
        // A voronoi diagram with just one reference point is just a loosely defined face.
        face := HEFace{
            ReferencePoint: points[0],
            EEdge:          EmptyEdge,
            IsClosed:       false,
        }
        return Voronoi{
            heVertexList: []HEVertex{},
            heEdgeList:   []HEEdge{},
            heFaceList:   []HEFace{face},
        }
    } else {
        // Simple voronoi with two reference points.
        p1 := points[0]
        p2 := points[1]
        bisector := PerpendicularBisector(p1, p2)

        edge1 := HEEdge {
            VOrigin: EmptyVertex,
            // ETwin empty for now
            ENext: EmptyEdge,
            // FFace empty for now
            TmpEdge: bisector,
        }
        edge2 := HEEdge {
            VOrigin: EmptyVertex,
            ETwin: &edge1,
            ENext: EmptyEdge,
            // FFace empty for now
            TmpEdge: bisector,
        }
        edge1.ETwin = &edge2


        face1 := HEFace {
            ReferencePoint: p1,
            EEdge: &edge1,
            IsClosed: false,
        }
        face2 := HEFace {
            ReferencePoint: p2,
            EEdge: &edge2,
            IsClosed: false,
        }

        // Set references of the faces.
        edge1.FFace = &face1
        edge2.FFace = &face2

        return Voronoi {
            heVertexList: []HEVertex{},
            heEdgeList:   []HEEdge{edge1, edge2},
            heFaceList:   []HEFace{face1, face2},
        }
    }

    fmt.Println("ERROR!")
    return Voronoi{}
}

// Calculates the face with the 'best' (y-coord) reference point.
// Can be used to calculate the lowest or highest point with the appropriate function.
func (ch ConvexHull) betterPoint(isBetter func(y1, y2 float32) bool) *HEFace {
    bestFace := ch[0]
    for _,face := range ch {
        if isBetter(face.ReferencePoint.Y, bestFace.ReferencePoint.Y) {
            bestFace = face
        }
    }
    return bestFace
}

// Calculates the highest intersection of the given bisector and any edge of the
// given face. With the restriction, that it can't be 'lastEdge'.
func calcHighestIntersection(bisector Edge, p *HEFace, lastEdge *HEEdge) (*HEEdge, Vector) {

    bisector.Amplify(100.0)

    // Try to go one direction. If we can't get all the way around, we go the other
    // direction as well!
    edge := p.EEdge
    line := edge.Line(true)
    bestIntersection := LineIntersection(bisector, line)
    bestEdge := edge
    edge = edge.ENext

    // Find be highest/best intersection of the bisector and and edge.
    for edge != EmptyEdge && edge != p.EEdge {

        line = edge.Line(true)
        intersection := LineIntersection(bisector, line)
        if intersection.Y > bestIntersection.Y {
            bestIntersection = intersection
            bestEdge = edge
        }
    }

    return bestEdge, bestIntersection
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func mergeVoronoi(left, right Voronoi) Voronoi {

    h1 := left.ConvexHull()
    h2 := right.ConvexHull()

    isHigher := func(y1, y2 float32) bool {
        return y1 > y2
    }
    p := h1.betterPoint(isHigher)
    q := h2.betterPoint(isHigher)

    isLower := func(y1, y2 float32) bool {
        return y1 < y2
    }
    h1Down := h1.betterPoint(isLower)
    h2Down := h2.betterPoint(isLower)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge

    bisector := PerpendicularBisector(p.ReferencePoint, q.ReferencePoint)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for *p != *h1Down && *q != *h2Down {

        edgeP, locationP := calcHighestIntersection(bisector, p, lastPEdge)
        edgeQ, locationQ := calcHighestIntersection(bisector, q, lastQEdge)

        // Which intersection point do we take?
        if locationP.Y > locationQ.Y {
            // Do something with the voronoi...
            p = edgeP.ETwin.FFace
        } else {
            // Do something with the voronoi...
            q = edgeQ.ETwin.FFace
        }
        bisector = PerpendicularBisector(p.ReferencePoint, q.ReferencePoint)
    }

    return Voronoi {
        heVertexList: append(left.heVertexList, right.heVertexList...),
        heEdgeList:   append(left.heEdgeList, right.heEdgeList...),
        heFaceList:   append(left.heFaceList, right.heFaceList...),
    }
}

// Voronoi divide and conquer entry point
func divideAndConquer(points PointList) Voronoi {

    l := len(points)
    // Recursion break on two points (what about 1 vs 3?)
    if l <= 2 {
        return createTrivialVoronoi(points)
    }

    left  := divideAndConquer(points[:l/2])
    right := divideAndConquer(points[l/2:])

    return mergeVoronoi(left, right)

}

// Sorts the points and returns a voronoi tessellation.
func CreateVoronoi(pointList PointList) Voronoi {
    sort.Sort(pointList)
    return divideAndConquer(pointList)
}

func draw(v Voronoi) {

}

func main() {
    var r = rand.New(rand.NewSource(0))
    var pointList PointList

    for i:=0; i < 10; i++ {
        pointList = append(pointList, Vector{r.Float32()*500., r.Float32()*500., 0})
    }

    var voronoi = CreateVoronoi(pointList)

    voronoi.pprint()
}




