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
                // Best case. We have both endpoints.
                case v1 != EmptyVertex && v2 != EmptyVertex:
                    edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                // We have the "left" endpoint.
                case v1 != EmptyVertex && v2 == EmptyVertex:
                    edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                // We have the "right" endpoint.
                case v1 == EmptyVertex && v2 != EmptyVertex:
                    edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                // We don't have any endpoints.
                case v1 == EmptyVertex && v2 == EmptyVertex:
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

    //bisector.Amplify(500.0)
    //bisector = Amplify(bisector, 50.0)

    lastV := lastVertex != EmptyVertex
    lastY := float32(0.0)
    if lastV {
        lastY = v.vertices[lastVertex].Pos.Y
    }

    // Try to go one direction. If we can't get all the way around, we go the other
    // direction as well!
    edge := v.faces[face].EEdge
    veryFirstEdge := edge

    // For a voronoi with only one or two faces and no edge...
    if edge == EmptyEdge {
        return EmptyEdge, Vector{}
    }

    firstIntersects, bestIntersection := LineIntersection4(bisector, createLine(v, edge, true))
    if firstIntersects && lastV && bestIntersection.Y > lastY {
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

            // For an intersection to be considered, it must satisfy the following conditions:
            // - Must intersect!
            // - Must be better, than the best one so far (if there is one)
            // - Must be below the last Vertex (if there is one!)
            if intersects && (!firstIntersects || intersection.Y > bestIntersection.Y) && (!lastV || (intersection.Y < lastY)) {

                firstIntersects = true
                bestIntersection = intersection
                bestEdge = edge
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

func commonSupportLine(v *Voronoi, h1, h2 ConvexHull, betterSide func(v1, v2 Vector) bool, betterPoint func(v1, v2, test Vector) bool) (FaceIndex, FaceIndex) {

    ai := h1.bestFace(v, betterSide)
    bi := h2.bestFace(v, betterSide)

    finished := false

    for !finished {
        finished = true
        // iterating through the convex hull points...
        for betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[(ai+1)%len(h1)]].ReferencePoint) {
            ai = (ai+1)%len(h1)
            finished = false
        }
        // iterating through the convex hull points...
        for betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[(bi+1)%len(h2)]].ReferencePoint) {
            bi = (bi+1)%len(h2)
            finished = false
        }
    }

    return h1[ai], h2[bi]
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func (v *Voronoi)mergeVoronoi(left, right VoronoiEntryFace) VoronoiEntryFace {

    voronoiEntry := left

    h1 := v.ConvexHull(left)
    h2 := v.ConvexHull(right)

    isRight := func(v1, v2 Vector) bool {
        return v1.X > v2.X
    }
    // Upper common support line!
    p, q := commonSupportLine(v, h1, h2, isRight, IsLeft2D)

    isLeft := func(v1, v2 Vector) bool {
        return v1.X < v2.X
    }
    // Lower common support line!
    h1Down, h2Down := commonSupportLine(v,  h1, h2, isLeft, IsRight2D)

    fmt.Printf("Upper common support line: %v, %v\n", v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
    fmt.Printf("Lower common support line: %v, %v\n", v.faces[h1Down].ReferencePoint, v.faces[h2Down].ReferencePoint)

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

    //fmt.Println("bisector: ", bisector)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge, lastVertex)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge, lastVertex)

        fmt.Printf("\nmerge -- last: %v, locP: %v, locQ: %v, edgeP: %v, edgeQ: %v \n", lastMerge, locationP, locationQ, edgeP, edgeQ)

        switch {

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case lastMerge || (edgeP == EmptyEdge && edgeQ == EmptyEdge):
                fmt.Println("--> primitiv or last merge.")

                otherWayBisector := bisector
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp   := v.createEdge(EmptyVertex, EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex,  heEdgeUp,  EmptyEdge,  q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                }

                v.faces[p].EEdge = heEdgeUp

                // For merging primitive voronois.
                if v.faces[q].EEdge == EmptyEdge {
                    v.faces[q].EEdge = heEdgeDown
                }
                //lastMerge = true

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && (locationP.Y >= locationQ.Y || edgeQ == EmptyEdge):
                fmt.Println("--> intersection with p")
                heVertex := v.createVertex(locationP)

                if v.edges[v.edges[edgeP].ETwin].VOrigin != EmptyVertex {
                    fmt.Println("============ FOUND THAT SHIT P")
                }

                if v.edges[edgeP].FFace != p {
                    fmt.Println("============ ASSERTION. WRONG P")
                }

                otherWayBisector := bisector
                //bisector.Dir = Mult(bisector.Dir, -1)
                otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

                fmt.Println("bisector: ", bisector)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, nextPEdge, p, otherWayBisector)
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, EmptyEdge, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown

                v.edges[edgeP].ENext = heEdgeUp
                v.edges[v.edges[edgeP].ETwin].VOrigin = heVertex

                // heEdgeDown could  now be the first edge of face q. But only, if this is the first cut.
                if lastDownEdge != EmptyEdge {
                    v.edges[lastDownEdge].ENext = heEdgeDown
                } else {
                    if v.faces[q].EEdge == EmptyEdge {
                        v.faces[q].EEdge = heEdgeDown
                    }
                }

                nextPEdge    = v.edges[edgeP].ETwin
                lastDownEdge = heEdgeDown
                lastVertex   = heVertex
                lastPEdge    = edgeP
                lastQEdge    = EmptyEdge

                p = v.edges[v.edges[edgeP].ETwin].FFace

                // Special case for situations where more than three edges meet.
                // This could relat in kind of invalid Delaunay triangulations. At least regarding a triangle.
                // This should now work for situations, where 4 edges meet. I have to re-examine for even more edges meeting...
                if edgeQ != EmptyEdge && Equal(locationP, locationQ) {
                    fmt.Println("=========== Special case. More than three edges meet.")
                    v.edges[heEdgeDown].ENext = edgeQ
                    v.edges[edgeQ].VOrigin = heVertex
                    lastDownEdge = v.edges[edgeQ].ETwin
                    q = v.edges[v.edges[edgeQ].ETwin].FFace
                }

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && ((locationQ.Y >= locationP.Y) || edgeP == EmptyEdge):
                fmt.Println("--> intersection with q")
                heVertex := v.createVertex(locationQ)

                if v.edges[edgeQ].VOrigin != EmptyVertex {
                    fmt.Println("============ FOUND THAT SHIT Q")
                }

                if v.edges[edgeQ].FFace != q {
                    fmt.Println("============ ASSERTION. WRONG Q")
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
                    voronoiEntry = VoronoiEntryFace(q)
                }

                // If p was a single face without edges
                if v.faces[p].EEdge == EmptyEdge {
                    v.faces[p].EEdge = heEdgeUp
                }

                v.edges[edgeQ].VOrigin = heVertex

                nextPEdge    = heEdgeUp
                lastDownEdge = v.edges[edgeQ].ETwin
                lastVertex   = heVertex
                lastQEdge    = edgeQ
                lastPEdge    = EmptyEdge

                q = v.edges[v.edges[edgeQ].ETwin].FFace
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
        bisector = Amplify(bisector, 50.0)

        if lastMerge {
            break
        }

    }


    fmt.Printf("FINISHED MERGE OF %v AND %v\n\n", left, right)

    return voronoiEntry
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
        pointList = append(pointList, Vector{r.Float32()*100., r.Float32()*100., 0})
    }

    // Weird shit happening
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

    // Works fine!

    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 30., 0})
    pointList = append(pointList, Vector{80., 20., 0})
    pointList = append(pointList, Vector{55., 40., 0})


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
    */

    // Looping infinitly. Problem occurs when the last point is added!
    /*
    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 30., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{25., 40., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    */

    v := CreateVoronoi(pointList)
    v.pprint()

    v.createImage()

}




