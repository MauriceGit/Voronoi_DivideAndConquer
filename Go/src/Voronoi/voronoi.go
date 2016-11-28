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

    bisector.Amplify(500.0)

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
    nextPEdge  := EmptyEdge
    lastDownEdge  := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := EmptyVertex

    bisector := PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
    bisector.Amplify(500.0)

    fmt.Println("bisector: ", bisector)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge, lastVertex)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge, lastVertex)

        fmt.Println("merge", lastMerge, locationP, locationQ, edgeP, edgeQ)

        switch {

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case lastMerge || (edgeP == EmptyEdge && edgeQ == EmptyEdge):
                fmt.Println("--> primitiv or last merge.")
                // nextPEdge is false. It should be something like: lastPUpEdge!!!
                heEdgeUp   := v.createEdge(EmptyVertex, EmptyEdge, nextPEdge, p, bisector)
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
                lastMerge = true

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && ((locationP.Y >= locationQ.Y) || edgeQ == EmptyEdge):
                fmt.Println("--> intersection with p")
                heVertex := v.createVertex(locationP, EmptyEdge)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, nextPEdge, p, bisector)

                // We don't know yet, what the next Edge is! So we have to set that later!
                heEdgeDown := v.createEdge(lastVertex, heEdgeUp, EmptyEdge, q, bisector)
                v.edges[heEdgeUp].ETwin = heEdgeDown
                v.vertices[heVertex].ELeaving = heEdgeUp

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

                p = v.edges[v.edges[edgeP].ETwin].FFace

                // Special case for situations where more than three edges meet.
                // This could relat in kind of invalid Delaunay triangulations. At least regarding a triangle.
                // This should now work for situations, where 4 edges meet. I have to re-examine for even more edges meeting...
                if edgeQ != EmptyEdge && Equal(locationP, locationQ) {
                    v.edges[heEdgeDown].ENext = edgeQ
                    v.edges[edgeQ].VOrigin = heVertex
                    lastDownEdge = v.edges[edgeQ].ETwin
                    q = v.edges[v.edges[edgeQ].ETwin].FFace
                }

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && ((locationQ.Y >= locationP.Y) || edgeP == EmptyEdge):
                fmt.Println("--> intersection with q")
                heVertex := v.createVertex(locationQ, EmptyEdge)

                heEdgeUp := v.createEdge(heVertex, EmptyEdge, nextPEdge, p, bisector)

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

                nextPEdge    = heEdgeUp
                lastDownEdge = v.edges[edgeQ].ETwin
                lastVertex   = heVertex
                lastQEdge    = edgeQ

                q = v.edges[v.edges[edgeQ].ETwin].FFace
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
        bisector.Amplify(500.0)

        if lastMerge {
            break
        }

        fmt.Println("bisector: ", bisector)
    }


    fmt.Printf("FINISHED MERGE OF %v AND %v\n", left, right)

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

type Circle struct {
    X, Y, R float64
}

func (c *Circle) Brightness(x, y float64) uint8 {
    var dx, dy float64 = c.X - x, c.Y - y
    d := math.Sqrt(dx*dx+dy*dy) / c.R
    if d > 1 {
        return 0
    } else {
        return 255
    }
}

func drawCircle(posX, posY, radius float64) {
    var w, h int = 1000, 1000
    hw, hh := posX, posY
    cr := &Circle{hw, hh, radius}
    //cg := &Circle{hw - radius*math.Sin(θ), hh - radius*math.Cos(θ), radius}
    //cb := &Circle{hw - radius*math.Sin(-θ), hh - radius*math.Cos(-θ), radius}

    m := image.NewRGBA(image.Rect(0, 0, w, h))
    for x := 0; x < w; x++ {
        for y := 0; y < h; y++ {
            c := color.RGBA{
                cr.Brightness(float64(x), float64(y)),
                cr.Brightness(float64(x), float64(y)),
                cr.Brightness(float64(x), float64(y)),
                255,
            }
            if cr.Brightness(float64(x), float64(y)) > 0 {
                m.Set(x, y, c)
            }
        }
    }

    f, err := os.OpenFile("rgb.png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

func draw(pointList PointList) {
    // Initialize the graphic context on an RGBA image
    dest := image.NewRGBA(image.Rect(0, 0, 500, 500.0))
    gc := draw2dimg.NewGraphicContext(dest)

    // Set some properties
    gc.SetFillColor(color.RGBA{0x44, 0xff, 0x44, 0xff})
    gc.SetStrokeColor(color.RGBA{0xff, 0x00, 0x00, 0xff})
    gc.SetLineWidth(5)

    for _,p := range pointList {

        p = p
        //gc.Set(2, 3, color.RGBA{255, 0, 0, 255})

    }

    // Save to file
    draw2dimg.SaveToPngFile("hello.png", dest)
}

func (v *Voronoi) draw() {
    // Initialize the graphic context on an RGBA image
    dest := image.NewRGBA(image.Rect(0, 0, 500, 500.0))
    gc := draw2dimg.NewGraphicContext(dest)

    // Set some properties
    gc.SetFillColor(color.RGBA{0x44, 0xff, 0x44, 0xff})
    gc.SetStrokeColor(color.RGBA{0xff, 0x00, 0x00, 0xff})
    gc.SetLineWidth(5)

    for i,_ := range v.edges {
        if v.edges[EdgeIndex(i)].VOrigin == EmptyVertex || v.edges[v.edges[EdgeIndex(i)].ETwin].VOrigin == EmptyVertex {
            continue
        }



        line := createLine(v, EdgeIndex(i), false)

        fmt.Println("Test? ", line)

        gc.MoveTo(float64(line.Pos.X*5.), 500-float64(line.Pos.Y*5.))
        gc.LineTo(float64(Add(line.Pos, line.Dir).X*5.), 500-float64(Add(line.Pos, line.Dir).Y*5.))
        gc.Close()
        gc.FillStroke()

    }

    gc.SetStrokeColor(color.RGBA{0xff, 0xff, 0x00, 0xff})

    for _,f := range v.faces {

        gc.MoveTo(float64(f.ReferencePoint.X*5.), 500-float64(f.ReferencePoint.Y*5.))
        gc.LineTo(float64(f.ReferencePoint.X*5.)+5, 500-float64(f.ReferencePoint.Y*5.)+5)
        gc.Close()
        gc.FillStroke()
    }

    gc.SetStrokeColor(color.RGBA{0x00, 0xff, 0x00, 0xff})

    for _,ve := range v.vertices {

        gc.MoveTo(float64(ve.Pos.X*5.), 500-float64(ve.Pos.Y*5.))
        gc.LineTo(float64(ve.Pos.X*5.)+5, 500-float64(ve.Pos.Y*5.)+5)
        gc.Close()
        gc.FillStroke()
    }

    // Save to file
    draw2dimg.SaveToPngFile("voronoi.png", dest)
}


func main() {

    var r = rand.New(rand.NewSource(0))
    var pointList PointList

    for i:= 0; i < 10; i++ {
        pointList = append(pointList, Vector{r.Float32()*100., r.Float32()*100., 0})
    }

    //draw(pointList)

    //v := CreateVoronoi(pointList)
    //v.pprint()
    //v.draw()

    drawCircle(200, 300, 50)

}




