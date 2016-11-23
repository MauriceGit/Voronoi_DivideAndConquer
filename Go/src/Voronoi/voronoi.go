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
    heVertexList    []*HEVertex
    heEdgeList      []*HEEdge
    heFaceList      []*HEFace
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
    var outerFace = EmptyFace

    // find outer face as starting point for convex hull.
    for i,face := range v.heFaceList {
        if !face.IsClosed {
            outerFace = v.heFaceList[i]
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
            for *lastEdge.ENext != *EmptyEdge {
                lastEdge = lastEdge.ENext
            }
            nextFace = lastEdge.ETwin.FFace
        }
    }

    return convexHullList
}

// So we can get a pointer of some data structure? What about scope issues?
func createFace(refPoint Vector, eEdge *HEEdge, closed bool) *HEFace {
    return &HEFace {
        ReferencePoint: refPoint,
        EEdge:          eEdge,
        IsClosed:       closed,
    }
}

// So we can get a pointer of some data structure? What about scope issues?
func createEdge(vOrigin *HEVertex, eTwin, eNext *HEEdge, fFace *HEFace, tmpEdge Edge) *HEEdge {
    return &HEEdge {
        VOrigin:    vOrigin,
        ETwin:      eTwin,
        ENext:      eNext,
        FFace:      fFace,
        TmpEdge:    tmpEdge,
    }
}

// So we can get a pointer of some data structure? What about scope issues?
func createVertex(pos Vector, eLeaving *HEEdge) *HEVertex {
    return &HEVertex {
        Pos:        pos,
        ELeaving:   eLeaving,
    }
}

// There can only be 1 or 2 points in the list!!! This HAS to be guaranteed!
func createTrivialVoronoi(points PointList) Voronoi {

    if len(points) == 1 {
        // A voronoi diagram with just one reference point is just a loosely defined face.
        return Voronoi{
            heVertexList: []*HEVertex{},
            heEdgeList:   []*HEEdge{},
            heFaceList:   []*HEFace{createFace(points[0], EmptyEdge, false)},
        }
    } else {
        // Simple voronoi with two reference points.
        p1 := points[0]
        p2 := points[1]
        bisector := PerpendicularBisector(p1, p2)

        edge1 := createEdge(EmptyVertex, EmptyEdge, EmptyEdge, EmptyFace, bisector)
        edge2 := createEdge(EmptyVertex, edge1, EmptyEdge, EmptyFace, bisector)
        edge1.ETwin = edge2

        face1 := createFace(p1, edge1, false)
        face2 := createFace(p2, edge2, false)

        // Set references of the faces.
        edge1.FFace = face1
        edge2.FFace = face2

        return Voronoi {
            heVertexList: []*HEVertex{},
            heEdgeList:   []*HEEdge{edge1, edge2},
            heFaceList:   []*HEFace{face1, face2},
        }
    }

    fmt.Println("ERROR!")
    return Voronoi{}
}

// Calculates the face with the 'best' (y-coord) reference point.
// Can be used to calculate the lowest or highest point with the appropriate function.
func (ch ConvexHull) bestPoint(isBetter func(y1, y2 float32) bool) *HEFace {
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
    for *edge != *EmptyEdge && *edge != *p.EEdge && *edge != *p.EEdge.ETwin {

        line = edge.Line(true)
        intersection := LineIntersection(bisector, line)
        if intersection.Y > bestIntersection.Y {
            bestIntersection = intersection
            bestEdge = edge
        }
    }

    // So
    //if *bestEdge.FFace == *p {
    //    bestEdge = bestEdge.ETwin
    //}

    return bestEdge, bestIntersection
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func mergeVoronoi(left, right Voronoi) Voronoi {

    mergedVoronoi := Voronoi {
        heVertexList: append(left.heVertexList, right.heVertexList...),
        heEdgeList:   append(left.heEdgeList, right.heEdgeList...),
        heFaceList:   append(left.heFaceList, right.heFaceList...),
    }

    h1 := left.ConvexHull()
    h2 := right.ConvexHull()

    isHigher := func(y1, y2 float32) bool {
        return y1 > y2
    }
    // p and q are faces!
    p := h1.bestPoint(isHigher)
    q := h2.bestPoint(isHigher)

    isLower := func(y1, y2 float32) bool {
        return y1 < y2
    }
    h1Down := h1.bestPoint(isLower)
    h2Down := h2.bestPoint(isLower)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge
    // This is either lastPEdge or lastQEdge.
    lastUpEdge  := EmptyEdge
    lastDownEdge  := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := EmptyVertex

    bisector := PerpendicularBisector(p.ReferencePoint, q.ReferencePoint)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for *p != *h1Down && *q != *h2Down {

        edgeP, locationP := calcHighestIntersection(bisector, p, lastPEdge)
        edgeQ, locationQ := calcHighestIntersection(bisector, q, lastQEdge)

        // Which intersection point do we take?
        if locationP.Y > locationQ.Y {

            heVertex := createVertex(locationP, EmptyEdge)

            heEdgeUp := createEdge(heVertex, EmptyEdge, lastUpEdge, p, bisector)

            // We don't know yet, what the next Edge is! So we have to set that later!
            heEdgeDown := createEdge(lastVertex, heEdgeUp, EmptyEdge, q, bisector)
            heEdgeUp.ETwin = heEdgeDown
            heVertex.ELeaving = heEdgeUp
            lastDownEdge.ENext = heEdgeDown

            // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
            if *lastDownEdge == *EmptyEdge {
                q.EEdge = heEdgeDown
            }

            mergedVoronoi.heVertexList = append(mergedVoronoi.heVertexList, heVertex)
            mergedVoronoi.heEdgeList = append(mergedVoronoi.heEdgeList, heEdgeUp)
            mergedVoronoi.heEdgeList = append(mergedVoronoi.heEdgeList, heEdgeDown)


            lastUpEdge   = heEdgeUp
            lastDownEdge = heEdgeDown
            lastVertex   = heVertex
            lastPEdge    = edgeP

            p = edgeP.ETwin.FFace

        } else {

            heVertex := createVertex(locationQ, EmptyEdge)

            // todo: the next edge is NOT lastUpEdge... Only sometimes... shit.
            heEdgeUp := createEdge(heVertex, EmptyEdge, lastUpEdge, p, bisector)

            // We don't know yet, what the next Edge is! So we have to set that later!
            // Shit. Next edge is also not the next edge every time...
            heEdgeDown := createEdge(lastVertex, heEdgeUp, EmptyEdge, q, bisector)
            heEdgeUp.ETwin = heEdgeDown
            heVertex.ELeaving = heEdgeUp
            lastDownEdge.ENext = heEdgeDown

            // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
            if *lastDownEdge == *EmptyEdge {
                q.EEdge = heEdgeDown
            }

            mergedVoronoi.heVertexList = append(mergedVoronoi.heVertexList, heVertex)
            mergedVoronoi.heEdgeList = append(mergedVoronoi.heEdgeList, heEdgeUp)
            mergedVoronoi.heEdgeList = append(mergedVoronoi.heEdgeList, heEdgeDown)


            lastUpEdge   = heEdgeUp
            lastDownEdge = heEdgeDown
            lastVertex   = heVertex
            lastPEdge    = edgeP

            p = edgeP.ETwin.FFace


            // Do something with the voronoi...
            q = edgeQ.ETwin.FFace
            //lastEdge = edgeQ
        }
        bisector = PerpendicularBisector(p.ReferencePoint, q.ReferencePoint)
    }

    // Combine the datasets. They should be connected now anyway.
    return mergedVoronoi
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




