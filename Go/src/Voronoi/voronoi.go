package main

import (
    . "Voronoi/Vector"
    . "Voronoi/HalfEdge"
    "fmt"
    "math/rand"
    "sort"
    //"math"
)

// HAS to be an infinite face!
type VoronoiEntryFace FaceIndex

// List of indides to the faces, that form the convex hull.
type ConvexHull []FaceIndex

type Voronoi struct {
    vertices           []HEVertex
    firstFreeVertexPos VertexIndex
    edges              []HEEdge
    firstFreeEdgePos   EdgeIndex
    faces              []HEFace
    firstFreeFacePos   FaceIndex
}


func (v *Voronoi) pprint () {
    fmt.Println("Voronoi:")
    fmt.Printf("   Vertices (%v): ", v.firstFreeVertexPos)
    var dummyV HEVertex
    for _,ve := range v.vertices {
        if ve != dummyV {
            fmt.Printf("%v, ", ve)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Edges    (%v): ", v.firstFreeEdgePos)
    var dummyE HEEdge
    for _,e := range v.edges {
        if e != dummyE {
            fmt.Printf("%v, ", e)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Faces    (%v): ", v.firstFreeFacePos)
    var dummyF HEFace
    for _,f := range v.faces {
        if f != dummyF {
            fmt.Printf("%v, ", f)
        }
    }
    fmt.Printf("\n")
}

// Calculates a convex hull for the given voronoi. Runs in O(n).
// The convex hull consists of a list of Reference points.
// As there is an exactly 1:1 representation of Reference points and HEFaces,
// a list of face indices is returned in Order (counter-clockwise) of the convex hull.
func (v *Voronoi)ConvexHull(fEntry VoronoiEntryFace) ConvexHull {

    f := FaceIndex(fEntry)
    convexHullList := []FaceIndex{FaceIndex(f)}

    edge := v.faces[f].EEdge

    // If the voronoi consists of just one refPoint, there is no edge.
    if edge != EmptyEdge {

        nextFace := v.edges[v.edges[edge].ETwin].FFace

        for nextFace != f {
            convexHullList = append(convexHullList, nextFace)

            nextFace = v.edges[v.edges[v.faces[nextFace].EEdge].ETwin].FFace
        }
    }

    return convexHullList
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createFace(refPoint Vector, eEdge EdgeIndex) FaceIndex {
    v.faces[v.firstFreeFacePos] = HEFace {
        ReferencePoint: refPoint,
        EEdge:          eEdge,
    }
    fmt.Println("face ", v.firstFreeFacePos, " == ", v.faces[v.firstFreeFacePos])
    v.firstFreeFacePos += 1
    return v.firstFreeFacePos-1
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createEdge(vOrigin VertexIndex, eTwin, eNext EdgeIndex, fFace FaceIndex, tmpEdge Edge) EdgeIndex {
    v.edges[v.firstFreeEdgePos] = HEEdge {
        VOrigin:    vOrigin,
        ETwin:      eTwin,
        ENext:      eNext,
        FFace:      fFace,
        TmpEdge:    tmpEdge,
    }
    v.firstFreeEdgePos += 1
    return v.firstFreeEdgePos-1
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createVertex(pos Vector, eLeaving EdgeIndex) VertexIndex {
    v.vertices[v.firstFreeVertexPos] = HEVertex {
        Pos:        pos,
        ELeaving:   eLeaving,
    }
    v.firstFreeVertexPos += 1
    return v.firstFreeVertexPos-1
}

// Calculates the face with the 'best' (y-coord) reference point.
// Can be used to calculate the lowest or highest point with the appropriate function.
func (ch ConvexHull) bestFace(v *Voronoi, isBetter func(y1, y2 float32) bool) FaceIndex {
    bestFace := ch[0]
    for _,face := range ch {
        if isBetter(v.faces[face].ReferencePoint.Y, v.faces[bestFace].ReferencePoint.Y) {
            bestFace = face
        }
    }
    return bestFace
}

// Creates a line from the HEEdge. Depeding on its state from the edge or the TmpEdge.
func createLine (v *Voronoi, e EdgeIndex, amplified bool) Edge {
    if v.edges[e].VOrigin == EmptyVertex || v.edges[v.edges[e].ETwin].VOrigin == EmptyVertex {
        if amplified {
            tmpE := v.edges[e].TmpEdge.Copy()
            tmpE.Amplify(100.0)
            return tmpE
        } else {
            return v.edges[e].TmpEdge
        }
    } else {
        return Edge {
            Pos: v.vertices[v.edges[e].VOrigin].Pos,
            Dir: Sub(v.vertices[v.edges[v.edges[e].ETwin].VOrigin].Pos, v.vertices[v.edges[e].VOrigin].Pos),
        }
    }
}

// Calculates the highest intersection of the given bisector and any edge of the
// given face. With the restriction, that it can't be 'lastEdge'.
func calcHighestIntersection(v *Voronoi, bisector Edge, face FaceIndex, lastEdge EdgeIndex) (EdgeIndex, Vector) {

    bisector.Amplify(100.0)

    // Try to go one direction. If we can't get all the way around, we go the other
    // direction as well!
    edge := v.faces[face].EEdge

    // For a voronoi with only one or two faces and no edge...
    if edge == EmptyEdge {
        return EmptyEdge, Vector{}
    }

    firstIntersects, bestIntersection := LineIntersection4(bisector, createLine(v, edge, true))
    bestEdge := edge
    edge = v.edges[edge].ENext

    // Find be highest/best intersection of the bisector and edge.
    for edge != EmptyEdge && edge != v.faces[face].EEdge && edge != v.edges[v.faces[face].EEdge].ETwin {

        intersects, intersection := LineIntersection4(bisector, createLine(v, edge, true))
        if intersects && (!firstIntersects || (intersection.Y > bestIntersection.Y)) {
            firstIntersects = true
            bestIntersection = intersection
            bestEdge = edge
        }
    }

    fmt.Println("The very best intersection iiiis: ", bestEdge, bestIntersection)

    return bestEdge, bestIntersection
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func (v *Voronoi)mergeVoronoi(left, right VoronoiEntryFace) VoronoiEntryFace {

    voronoiEntry := left

    h1 := v.ConvexHull(left)
    h2 := v.ConvexHull(right)

    isHigher := func(y1, y2 float32) bool {
        return y1 > y2
    }
    // p and q are faces!
    p := h1.bestFace(v, isHigher)
    q := h2.bestFace(v, isHigher)

    isLower := func(y1, y2 float32) bool {
        return y1 < y2
    }
    h1Down := h1.bestFace(v, isLower)
    h2Down := h2.bestFace(v, isLower)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge
    // This is either lastPEdge or lastQEdge.
    lastUpEdge  := EmptyEdge
    lastDownEdge  := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := EmptyVertex

    bisector := PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge)

        fmt.Println("merge", lastMerge, locationP, locationQ, edgeP, edgeQ)

        switch {

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case lastMerge || (edgeP == EmptyEdge && edgeQ == EmptyEdge):
                // lastUpEdge is false. It should be something like: lastPUpEdge!!!
                heEdgeUp   := v.createEdge(EmptyVertex, EmptyEdge, lastUpEdge, p, bisector)
                heEdgeDown := v.createEdge(lastVertex,  heEdgeUp,  EmptyEdge,  q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                v.faces[p].EEdge = heEdgeUp

                // For merging primitive voronois.
                if edgeP == EmptyEdge && edgeQ == EmptyEdge {
                    v.faces[q].EEdge = heEdgeDown
                }

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && ((locationP.Y > locationQ.Y) || edgeQ == EmptyEdge):
                heVertex := v.createVertex(locationP, EmptyEdge)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, lastUpEdge, p, bisector)

                // We don't know yet, what the next Edge is! So we have to set that later!
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, EmptyEdge, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown
                v.vertices[heVertex].ELeaving = heEdgeUp
                v.edges[lastDownEdge].ENext = heEdgeDown

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                }

                // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                if lastDownEdge == EmptyEdge {
                    v.faces[q].EEdge = heEdgeDown
                }

                lastUpEdge   = heEdgeUp
                lastDownEdge = heEdgeDown
                lastVertex   = heVertex
                lastPEdge    = edgeP

                p = v.edges[v.edges[edgeP].ETwin].FFace

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && ((locationQ.Y > locationP.Y) || edgeP == EmptyEdge):
                heVertex := v.createVertex(locationQ, EmptyEdge)

                // todo: the next edge is NOT lastUpEdge... Only sometimes... shit.
                heEdgeUp := v.createEdge(heVertex, EmptyEdge, lastUpEdge, p, bisector)

                // Next edge is always edgeQ!!! Has to be!
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, edgeQ, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown
                v.vertices[heVertex].ELeaving = heEdgeUp

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                } else {
                    // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                    v.faces[q].EEdge = heEdgeDown
                    voronoiEntry = VoronoiEntryFace(q)
                }

                // If p was a single face without edges
                if v.faces[p].EEdge == EmptyEdge {
                    v.faces[p].EEdge = heEdgeUp
                }

                v.edges[edgeQ].VOrigin = heVertex

                v.edges[heEdgeDown].ENext = edgeQ


                lastUpEdge   = heEdgeUp
                lastDownEdge = v.edges[edgeQ].ETwin
                lastVertex   = heVertex
                lastQEdge    = edgeQ

                q = v.edges[v.edges[edgeQ].ETwin].FFace
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)

        if lastMerge {
            break
        }
    }

    return voronoiEntry
}

// Voronoi divide and conquer entry point
func (v *Voronoi)divideAndConquer(points PointList) VoronoiEntryFace {
    l := len(points)

    if l == 1 {
        // This is, by definition, an outer face.
        return VoronoiEntryFace(v.createFace(points[0], EmptyEdge))
    }

    left  := v.divideAndConquer(points[:l/2])
    right := v.divideAndConquer(points[l/2:])

    return v.mergeVoronoi(left, right)
}

// Sorts the points and returns a voronoi tessellation.
func CreateVoronoi(pointList PointList) Voronoi {
    sort.Sort(pointList)
    n := len(pointList)

    // See: http://www.cs.wustl.edu/~pless/546/lectures/L11.html
    // for the calculations of maximum voronoi object count.
    v := Voronoi {
        vertices:           make([]HEVertex, 2*n-5 + 3),
        firstFreeVertexPos: 0,
        edges:              make([]HEEdge, 2*(3*n-6) + 6),
        firstFreeEdgePos:   0,
        faces:              make([]HEFace, n),
        firstFreeFacePos:   0,
    }

    v.divideAndConquer(pointList)

    return v
}

func (v Voronoi)draw() {

}

func main() {
    var r = rand.New(rand.NewSource(0))
    var pointList PointList

    for i:=0; i < 0; i++ {
        pointList = append(pointList, Vector{r.Float32()*500., r.Float32()*500., 0})
    }

    pointList = append(pointList, Vector{0,0,0})
    pointList = append(pointList, Vector{10,0,0})

    pointList = append(pointList, Vector{5,5,0})

    v := CreateVoronoi(pointList)
    v.pprint()

}




