package main

import (
    . "Voronoi/Vector"
    . "Voronoi/HalfEdge"
    "fmt"
    "math"
    "os"
    "math/rand"
    "sort"
    //"math"
    "github.com/llgcode/draw2d/draw2dimg"
    "image"
    "image/color"
    "image/png"
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

////////////////////////////////////////////////////////////////////////
//  Pretty Print the  Voronoi Attributes
////////////////////////////////////////////////////////////////////////

func (v *Voronoi) pprint () {
    fmt.Println("Voronoi:")
    fmt.Printf("   Vertices (%v):\n", v.firstFreeVertexPos)
    var dummyV HEVertex
    for i,ve := range v.vertices {
        if ve != dummyV {
            fmt.Printf("\t%v:\tPos: %v\n", i, ve.Pos)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Edges    (%v):\n", v.firstFreeEdgePos)
    var dummyE HEEdge
    for i,e := range v.edges {
        if e != dummyE {
            fmt.Printf("\t%v:\tOrigin: %v,\tTwin: %v,\tNext: %v,\tFace: %v\n", i, e.VOrigin, e.ETwin, e.ENext, e.FFace)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Faces    (%v):\n", v.firstFreeFacePos)
    var dummyF HEFace
    for i,f := range v.faces {
        if f != dummyF {
            fmt.Printf("\t%v:\tRefPoint: %v,\tEdge:%v\n", i, f.ReferencePoint, f.EEdge)
        }
    }
    fmt.Printf("\n")
}

////////////////////////////////////////////////////////////////////////
//  Draw an Image of the Voronoi Tessellation
////////////////////////////////////////////////////////////////////////

type Circle struct {
    X, Y, R float64
}

func (c *Circle) insideCircle(x, y float64) bool {
    var dx, dy float64 = c.X - x, c.Y - y
    d := math.Sqrt(dx*dx+dy*dy) / c.R
    return d <= 1
}

func drawCircle(m *image.RGBA, posX, posY, radius int, c color.RGBA) {
    cr := &Circle{float64(posX), float64(posY), float64(radius)}

    for x := posX-radius; x < posX+radius; x++ {
        for y := posY-radius; y < posY+radius; y++ {
            if cr.insideCircle(float64(x), float64(y)) {
                m.Set(x, y, c)
            }
        }
    }
}

func (v *Voronoi)createImage() {
    var w, h int = 1000, 1000

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // Edges
    c := color.RGBA{0,0,255,255}
    gc := draw2dimg.NewGraphicContext(m)
    gc.SetStrokeColor(c)
    gc.SetLineWidth(3)
    for i,e := range v.edges {
        var tmp HEEdge
        if e != tmp {
            e2 := v.edges[e.ETwin]
            edge := Edge{}
            v1 := e.VOrigin
            v2 := e2.VOrigin
            switch {
                // Best case. We have both endpoints. And none of them Infinity ones.
                case v1 != EmptyVertex && v2 != EmptyVertex && v.vertices[v1].Pos != InfinitePoint && v.vertices[v2].Pos != InfinitePoint:
                    edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                // We have the "left" endpoint.
                case v1 != EmptyVertex && v.vertices[v1].Pos != InfinitePoint && v2 == EmptyVertex:
                    edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                // We have the "right" endpoint.
                case v1 == EmptyVertex && v2 != EmptyVertex && v.vertices[v2].Pos != InfinitePoint:
                    edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                // We don't have any endpoints.
                default:
                    // amplified line. Exceeding all boundaries. Infinite line.
                    edge = createLine(v, EdgeIndex(i), false)
            }

            gc.MoveTo(float64(edge.Pos.X*10), float64(h)-float64(edge.Pos.Y*10))
            gc.LineTo(float64(Add(edge.Pos, edge.Dir).X*10), float64(h)-float64(Add(edge.Pos, edge.Dir).Y*10))
            gc.FillStroke()
            gc.Close()
        }
    }

    // Faces/Reference Points!
    c = color.RGBA{255,0,0,255}
    for _,f := range v.faces {
        drawCircle(m, int(f.ReferencePoint.X*10), h-int(f.ReferencePoint.Y*10), 5, c)
    }

    // Vertices between edges
    c = color.RGBA{0,255,0,255}
    for _,ve := range v.vertices {
        var tmp HEVertex
        if ve != tmp {
            drawCircle(m, int(ve.Pos.X*10), h-int(ve.Pos.Y*10), 5, c)
        }
    }

    f, err := os.OpenFile("voronoi.png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

////////////////////////////////////////////////////////////////////////
//  Calculate a Voronoi Tessellation
////////////////////////////////////////////////////////////////////////


// Calculates a convex hull for the given voronoi. Runs in O(n).
// The convex hull consists of a list of Reference points.
// As there is an exactly 1:1 representation of Reference points and HEFaces,
// a list of face indices is returned in Order (counter-clockwise) of the convex hull.
func (v *Voronoi)ConvexHull(fEntry VoronoiEntryFace) ConvexHull {
    fmt.Println("Convex Hull:")
    startFace := FaceIndex(fEntry)
    convexHullList := []FaceIndex{FaceIndex(startFace)}

    finished := false

    // If the voronoi consists of just one refPoint, there is no edge.
    if v.faces[startFace].EEdge != EmptyEdge {
        lastFace := startFace
        nextFace := v.edges[v.edges[v.faces[startFace].EEdge].ETwin].FFace
        //nextEdge := v.edges[v.faces[nextFace].EEdge].ENext
        convexHullList = append(convexHullList, nextFace)

        for {
            lastFace = nextFace
            nextFace = v.edges[v.edges[v.faces[nextFace].EEdge].ETwin].FFace

            tmpE := v.edges[v.faces[lastFace].EEdge].ENext
            fmt.Printf("nextFace: %v\n", nextFace)
            switch {
                case nextFace == startFace && tmpE != EmptyEdge && v.edges[v.edges[tmpE].ETwin].FFace != startFace && lastFace == startFace:
                    // Normally, this would go back to the last face.
                    // Except, when we have a linear voronoi! If this doesn't go back, we continue this way!!!
                    fmt.Println("YES")
                    nextFace = v.edges[v.edges[tmpE].ETwin].FFace
                case nextFace != startFace && tmpE == EmptyEdge:
                    fmt.Println("this shit")
                    // End of linear voronoi. We don't get any further...
                    finished = true
                case nextFace == startFace:
                    fmt.Println("next == start")
                    // We are finished the normal way.
                    finished = true

            }

            if finished {
                break
            }

            convexHullList = append(convexHullList, nextFace)
        }
    }
    fmt.Println("Convex Hull: Done")
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
func (v *Voronoi)createVertex(pos Vector) VertexIndex {
    v.vertices[v.firstFreeVertexPos] = HEVertex {
        Pos:        pos,
    }
    v.firstFreeVertexPos += 1
    return v.firstFreeVertexPos-1
}

// Calculates the face with the 'best' (y-coord) reference point.
// Can be used to calculate the lowest or highest point with the appropriate function.
func (ch ConvexHull) bestFace(v *Voronoi, isBetter func(v1, v2 Vector) bool) int {
    bestFace := 0
    for i,face := range ch {
        if isBetter(v.faces[face].ReferencePoint, v.faces[ch[bestFace]].ReferencePoint) {
            bestFace = i
        }
    }
    return bestFace
}

// Creates a line from the HEEdge. Depeding on its state from the edge or the TmpEdge.
func createLine (v *Voronoi, e EdgeIndex, amplified bool) Edge {
    if v.edges[e].VOrigin == EmptyVertex || v.edges[v.edges[e].ETwin].VOrigin == EmptyVertex {
        if amplified {
            tmpE := v.edges[e].TmpEdge.Copy()
            tmpE.Amplify(500.0)
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
func calcHighestIntersection(v *Voronoi, bisector Edge, face FaceIndex, lastEdge EdgeIndex, lastVertex VertexIndex) (EdgeIndex, Vector) {

    lastV := lastVertex != EmptyVertex
    lastY := 0.0
    if lastV {
        lastY = v.vertices[lastVertex].Pos.Y
    }

    edge := v.faces[face].EEdge
    veryFirstEdge := edge

    // For a voronoi with only one face and no edge...
    if edge == EmptyEdge {
        return EmptyEdge, Vector{}
    }

    firstIntersects, bestIntersection := LineIntersection4(bisector, createLine(v, edge, true))
    if firstIntersects && lastV && bestIntersection != InfinitePoint && bestIntersection.Y > lastY {
        firstIntersects = false
    }

    if lastVertex != EmptyVertex && firstIntersects && Equal(bestIntersection, v.vertices[lastVertex].Pos) {
        firstIntersects = false
    }

    bestEdge := edge
    if !firstIntersects {
        bestEdge = EmptyEdge
    }
    if edge == lastEdge || lastEdge != EmptyEdge && edge == v.edges[lastEdge].ETwin {
        firstIntersects = false
        bestEdge = EmptyEdge
    }
    edge = v.edges[edge].ENext

    // Find be highest/best intersection of the bisector and edge.
    for edge != EmptyEdge && edge != veryFirstEdge {

        if edge != lastEdge && (lastEdge == EmptyEdge || edge != v.edges[lastEdge].ETwin) {
            intersects, intersection := LineIntersection4(bisector, createLine(v, edge, true))

            fmt.Printf("equal: %v, %v --> %v\n", intersection, v.vertices[lastVertex].Pos, Equal(intersection, v.vertices[lastVertex].Pos))

            // For an intersection to be considered, it must satisfy the following conditions:
            // - Must intersect!
            // - Must be better, than the best one so far (if there is one)
            // - Must be below the lastVertex (if there is one!)
            if intersects && (lastVertex == EmptyVertex || !Equal(intersection, v.vertices[lastVertex].Pos)) {
                // lower than the last vertex if there is one!
                if (!lastV || (intersection.Y < lastY-EPS)) &&
                    // higher than the last intersection, if there is one of any one, if the best one is Infinity.
                    (!firstIntersects || bestIntersection == InfinitePoint || (intersection != InfinitePoint && intersection.Y > bestIntersection.Y)) {

                    firstIntersects = true
                    bestIntersection = intersection
                    bestEdge = edge
                }
            }
        }
        edge = v.edges[edge].ENext
    }

    return bestEdge, bestIntersection
}

func upperCommonSupportLine(v *Voronoi, h1, h2 ConvexHull) (FaceIndex, FaceIndex) {
    isRight := func(v1, v2 Vector) bool {
        return v1.X > v2.X
    }

    ai := h1.bestFace(v, isRight)
    bi := h2.bestFace(v, isRight)

    finished := false

    for !finished {
        finished = true
        // iterating through the convex hull points...
        for IsLeft2D(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[(ai+1)%len(h1)]].ReferencePoint) {
            ai = (ai+1)%len(h1)
            finished = false
        }
        // iterating through the convex hull points...
        for IsLeft2D(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[(bi+1)%len(h2)]].ReferencePoint) {
            bi = (bi+1)%len(h2)
            finished = false
        }
    }

    return h1[ai], h2[bi]

}

// Calculates the upper or lower common support line for the two given convex hulls.
// Depending on the betterSide and betterPoint-function, it could potentially also calculate
// left/right common support lines, depending on the divide-and-conquer approach of the main algorithm.
func commonSupportLine(v *Voronoi, h1, h2 ConvexHull, betterSide func(v1, v2 Vector) bool, betterPoint func(v1, v2, test Vector) bool) (FaceIndex, FaceIndex) {

    ai := h1.bestFace(v, betterSide)
    bi := h2.bestFace(v, betterSide)

    finished := false

    for !finished {
        finished = true
        // iterating through the convex hull points...
        for len(h1) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[(ai+1)%len(h1)]].ReferencePoint) {
            ai = (ai+1)%len(h1)
            finished = false
        }
        // iterating through the convex hull points...
        for len(h2) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[(bi+1)%len(h2)]].ReferencePoint) {
            bi = (bi+1)%len(h2)
            finished = false
        }
    }

    return h1[ai], h2[bi]
}
// Helperfunctions for calculating the common support line.
// Were anonymous functions, but sourced out because they where getting too long.
func isBetterUp(v1, v2, test Vector) bool {
    if IsLeft2D(v1, v2, test) {
        return true
    }
    side := SideOfLine(v1, v2, test)
    if side <= EPS && side >= -EPS {
        fmt.Println("ON LINE!")
        return test.X >= v1.X && test.X <= v2.X
    }
    return false
}
func isBetterDown(v1, v2, test Vector) bool {
    if IsRight2D(v1, v2, test) {
        return true
    }
    side := SideOfLine(v1, v2, test)
    if side <= EPS && side >= -EPS {
        fmt.Println("ON LINE!")
        return test.X >= v1.X && test.X <= v2.X
    }
    return false
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func (v *Voronoi)mergeVoronoi(left, right VoronoiEntryFace) VoronoiEntryFace {

    h1 := v.ConvexHull(left)
    h2 := v.ConvexHull(right)

    isRight := func(v1, v2 Vector) bool {
        return v1.X > v2.X
    }
    // Upper common support line!
    p, q := commonSupportLine(v, h1, h2, isRight, isBetterUp)

    isLeft := func(v1, v2 Vector) bool {
        return v1.X <= v2.X
    }
    // Lower common support line!
    h1Down, h2Down := commonSupportLine(v,  h1, h2, isLeft, isBetterDown)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge
    // This is either lastPEdge or lastQEdge.
    nextPEdge  := EmptyEdge
    lastDownEdge  := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := EmptyVertex

    bisector := PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
    bisector = Amplify(bisector, 50.0)
    //bisector.Amplify(500.0)

    fmt.Printf("h1Down: %v, h2Down: %v\n", v.faces[h1Down].ReferencePoint, v.faces[h2Down].ReferencePoint)
    fmt.Printf("CH-Left: %v, CH-Right: %v (left: %v, right: %v)\n", h1, h2, left, right)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge, lastVertex)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge, lastVertex)

        fmt.Printf("\nmerge -- last: %v, locP: %v, locQ: %v, edgeP: %v, edgeQ: %v, lastVertex: %v \n", lastMerge, locationP, locationQ, edgeP, edgeQ, lastVertex)

        switch {

            // Infinite vertex Q
            case lastMerge && edgeQ != EmptyEdge && locationQ == InfinitePoint:
                fmt.Println("found Infinity at Q")

                heVertex := v.createVertex(locationQ)

                otherWayBisector := bisector
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp   := v.createEdge(heVertex,    EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex,  heEdgeUp,  edgeQ,     q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                v.edges[edgeQ].VOrigin = heVertex

                if v.faces[p].EEdge == EmptyEdge {
                    v.faces[p].EEdge = heEdgeUp
                }

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                } else {
                    // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                    v.faces[q].EEdge = heEdgeDown
                }

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case edgeP == EmptyEdge && edgeQ == EmptyEdge:
            // Last Merge or we merge primitive voronois or there is no intersection any more (shouldn't and can't happen!!!)
            //case lastMerge:

                fmt.Println("last merge")

                if edgeP != EmptyEdge && locationP == InfinitePoint {
                    fmt.Println("found it P")
                }

                otherWayBisector := bisector
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp   := v.createEdge(EmptyVertex, EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex,  heEdgeUp,  edgeQ,     q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                }

                // For merging primitive voronois.
                v.faces[p].EEdge = heEdgeUp

                if v.faces[q].EEdge == EmptyEdge {
                    v.faces[q].EEdge = heEdgeDown
                }

            // Equal Intersection with both edges. So 4 edges meet.
            case edgeQ != EmptyEdge && edgeP != EmptyEdge && Equal(locationP, locationQ) && (lastVertex == EmptyVertex || !Equal(locationP, v.vertices[lastVertex].Pos)):
                // This could result in kind of an invalid Delaunay triangulations. At least regarding a triangle.
                // This should now work for situations, where 4 edges meet. I have to re-examine for even more edges meeting...
                fmt.Println("=========== Special case. More than three edges meet.")

                heVertex := v.createVertex(locationP)
                otherWayBisector := bisector
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                heEdgeUp   := v.createEdge(heVertex,   EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp,  edgeQ,     q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                v.edges[edgeP].ENext = heEdgeUp
                v.edges[v.edges[edgeP].ETwin].VOrigin = heVertex

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                }

                v.edges[edgeQ].VOrigin = heVertex


                nextPEdge    = v.edges[edgeP].ETwin
                lastDownEdge = v.edges[edgeQ].ETwin
                lastVertex   = heVertex
                lastPEdge    = edgeP
                lastQEdge    = edgeQ

                p = v.edges[v.edges[edgeP].ETwin].FFace
                q = v.edges[v.edges[edgeQ].ETwin].FFace

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && (locationP.Y >= locationQ.Y || edgeQ == EmptyEdge) && (lastVertex == EmptyVertex || !Equal(locationP, v.vertices[lastVertex].Pos)):
                fmt.Println("--> intersection with p")
                heVertex := v.createVertex(locationP)

                vertex := v.edges[v.edges[edgeP].ETwin].VOrigin
                if vertex != EmptyVertex {
                    fmt.Println("----> Found a vertex that has do be deleted (P)")

                    var emptyV HEVertex

                    // Otherwise, we already deleted that vertex and
                    // the VOrigin is going to be overwritten anyway.
                    if v.vertices[vertex] != emptyV {
                        var emptyE HEEdge
                        // If we have an infinity-point, there is no edge after that.
                        if v.vertices[vertex].Pos != InfinitePoint {
                            // Delete his twin of the following edge of edgeP
                            v.edges[v.edges[v.edges[edgeP].ENext].ETwin] = emptyE
                            // Delete the following edge of edgeP
                            v.edges[v.edges[edgeP].ENext] = emptyE
                        }
                        // Delete the vertex.
                        v.vertices[v.edges[v.edges[edgeP].ETwin].VOrigin] = emptyV
                    }

                }

                otherWayBisector := bisector
                //bisector.Dir = Mult(bisector.Dir, -1)
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, nextPEdge, p, otherWayBisector)
                // IMPORTANT: by setting v.faces[q].EEdge as the following edge, we create a temporary (wrong) next-edge.
                // This will be overwritten in the next loop but with this, we don't lose a possible reference to the other edge of q.
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, v.faces[q].EEdge, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                v.edges[edgeP].ENext = heEdgeUp
                v.edges[v.edges[edgeP].ETwin].VOrigin = heVertex

                // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                } else {
                    //if v.faces[q].EEdge == EmptyEdge {
                        v.faces[q].EEdge = heEdgeDown
                    //}
                }

                if v.edges[edgeP].VOrigin == EmptyVertex {
                    v.faces[p].EEdge = edgeP
                }

                nextPEdge    = v.edges[edgeP].ETwin
                lastDownEdge = heEdgeDown
                lastVertex   = heVertex
                lastPEdge    = edgeP
                lastQEdge    = EmptyEdge

                p = v.edges[v.edges[edgeP].ETwin].FFace

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && ((locationQ.Y >= locationP.Y) || edgeP == EmptyEdge) && (lastVertex == EmptyVertex || !Equal(locationQ, v.vertices[lastVertex].Pos)):
                fmt.Println("--> intersection with q")
                heVertex := v.createVertex(locationQ)

                vertex := v.edges[edgeQ].VOrigin
                if vertex != EmptyVertex {
                    fmt.Println("----> Found a vertex that has do be deleted (Q) - ", v.vertices[vertex].Pos)


                    var emptyV HEVertex

                    // Otherwise, we already deleted that vertex and
                    // the VOrigin is going to be overwritten anyway.
                    if v.vertices[vertex] != emptyV {
                        var emptyE HEEdge

                        // If we have an infinity-point, there is no edge after that.
                        if v.vertices[vertex].Pos != InfinitePoint {
                            // We have to get the edge after that vertex and its twin. A bit tricky here...
                            v.edges[v.edges[v.edges[v.edges[v.edges[v.edges[edgeQ].ETwin].ENext].ETwin].ENext].ETwin] = emptyE

                            // Delete the other one
                            v.edges[v.edges[v.edges[v.edges[v.edges[edgeQ].ETwin].ENext].ETwin].ENext] = emptyE
                        }
                        // Delete the vertex.
                        v.vertices[v.edges[edgeQ].VOrigin] = emptyV
                    }
                }

                otherWayBisector := bisector
                //bisector.Dir = Mult(bisector.Dir, -1)
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, edgeQ, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                } else {
                    // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                    v.faces[q].EEdge = heEdgeDown
                    //voronoiEntry = VoronoiEntryFace(q)
                }

                // If p was a single face without edges
                if v.faces[p].EEdge == EmptyEdge {
                    //v.faces[p].EEdge = heEdgeUp
                }

                v.edges[edgeQ].VOrigin = heVertex

                nextPEdge    = heEdgeUp
                lastDownEdge = v.edges[edgeQ].ETwin
                lastVertex   = heVertex
                lastQEdge    = edgeQ
                lastPEdge    = EmptyEdge

                q = v.edges[v.edges[edgeQ].ETwin].FFace

        }

        if lastMerge {
            break
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
        bisector = Amplify(bisector, 50.0)

    }

    fmt.Printf("FINISHED MERGE OF %v AND %v\n\n", left, right)

    return left
}

// Voronoi divide and conquer entry point
func (v *Voronoi)divideAndConquer(points PointList) VoronoiEntryFace {
    l := len(points)

    // Recursion stops at one point!
    if l == 1 {
        // One face without anything else is the most primitive, trivial Voronoi tessellation!
        return VoronoiEntryFace(v.createFace(points[0], EmptyEdge))
    }

    return v.mergeVoronoi(v.divideAndConquer(points[:l/2]), v.divideAndConquer(points[l/2:]))
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

func main() {

    var r = rand.New(rand.NewSource(0))
    var pointList PointList

    for i:= 0; i < 0; i++ {
        pointList = append(pointList, Vector{r.Float64()*100., r.Float64()*100., 0})
    }

    // Weird shit happening
    // Endless loop with the last point removed...
    /*
    pointList = append(pointList, Vector{10., 15., 0})
    pointList = append(pointList, Vector{15., 65., 0})
    pointList = append(pointList, Vector{20., 40., 0})
    pointList = append(pointList, Vector{40., 5., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 60., 0})
    pointList = append(pointList, Vector{70., 70., 0})
    pointList = append(pointList, Vector{90., 10., 0})
    */

    // Special case. Intersection with both p and q equally.
    // Requires special handling!
    /*
    pointList = append(pointList, Vector{28, 30, 0})
    pointList = append(pointList, Vector{30, 20, 0})
    pointList = append(pointList, Vector{30, 40, 0})

    pointList = append(pointList, Vector{42, 30, 0})
    pointList = append(pointList, Vector{40, 20, 0})
    pointList = append(pointList, Vector{40, 40, 0})
    */

    // Linear dependent minimal points
    pointList = append(pointList, Vector{20, 10, 0})
    pointList = append(pointList, Vector{30, 10, 0})
    pointList = append(pointList, Vector{40, 10, 0})
    pointList = append(pointList, Vector{50, 10, 0})


    // Works fine!
    /*
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 30., 0})
    pointList = append(pointList, Vector{80., 20., 0})
    pointList = append(pointList, Vector{55., 40., 0})
    */

    // Wrong vertices and not sure about how to calculate the actual tessellation myself!
    // HA, Works now!
    /*
    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 50., 0})
    */

    // Wrong vertex. Minimum configuration that fails.
    // Update: Works now!
    /*pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{45., 45., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    */

    // Linear dependent points. This is currently a problem!
    /*
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 10., 0})
    */

    // Looping infinitly. Problem occurs when the last point is added!

    //pointList = append(pointList, Vector{10., 10., 0})
    //pointList = append(pointList, Vector{20., 20., 0})

    //pointList = append(pointList, Vector{25., 40., 0})

    //pointList = append(pointList, Vector{30., 10., 0})

    //pointList = append(pointList, Vector{40., 30., 0})
    //pointList = append(pointList, Vector{50., 20., 0})
    //pointList = append(pointList, Vector{60., 10., 0})

    // Special case where the intersection is equal to p and q.
    // Works now for this case!
    /*
    pointList = append(pointList, Vector{30., 50., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{60., 20., 0})
    pointList = append(pointList, Vector{70., 50., 0})
    */

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage()

}




