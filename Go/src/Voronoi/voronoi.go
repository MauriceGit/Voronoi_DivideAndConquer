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
    "time"
    "strconv"
    "errors"
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

type ChainElem struct {
    intersection    Vector
    edgeP           EdgeIndex
    edgeQ           EdgeIndex
    p               FaceIndex
    q               FaceIndex
    bisector        Edge
}


var g_recursions int = 0
var g_edgeIndexList []EdgeIndex

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

func (v *Voronoi)createImage(filename string, whole bool) {
    var w, h int = 1000, 1000

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // Edges
    // yellow
    //c := color.RGBA{255,255,0,255}
    // blue
    c := color.RGBA{0,0,255,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(2)
    if whole {
        gc.SetLineWidth(2)
        c = color.RGBA{255,255,0,255}
    }
    gc.SetStrokeColor(c)

    for i,e := range v.edges {
        var tmp HEEdge
        e2   := v.edges[e.ETwin]
        if e != tmp && e2 != tmp {

            edge := Edge{}
            v1   := e.VOrigin
            v2   := e2.VOrigin
            switch {
                // Best case. We have both endpoints. And none of them Infinity ones.
                case v1.Valid() && v2.Valid() && v.vertices[v1].Pos != InfinitePoint && v.vertices[v2].Pos != InfinitePoint:
                    //fmt.Printf("Edge in question: \n\t%v\n\t%v\n\t%v, %v\n", e, e2, v1, v2)
                    edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                // We have the "left" endpoint.
                case v1.Valid() && v.vertices[v1].Pos != InfinitePoint && !v2.Valid():
                    edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                // We have the "right" endpoint.
                case !v1.Valid() && v2.Valid() && v.vertices[v2].Pos != InfinitePoint:
                    edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                // We don't have any endpoints.
                default:
                    //i = i
                    // amplified line. Exceeding all boundaries. Infinite line.
                    edge = createLine(v, EdgeIndex(i), false)
            }



            if edge == (Edge{}) {
                //continue
            }

            if i == 10 || i == 11 {
                fmt.Printf("coloring now for something\n")
                tmpC := color.RGBA{0,0,0,255}
                gc.SetStrokeColor(tmpC)
            } else {
                gc.SetStrokeColor(c)
            }

            //fmt.Printf("From %v|%v ----> %v|%v\n", edge.Pos.X*10., float64(h) - edge.Pos.Y*10., Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)

            gc.MoveTo(edge.Pos.X*10., float64(h) - edge.Pos.Y*10.)
            gc.LineTo(Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)
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
            drawCircle(m, int(ve.Pos.X*10), h-int(ve.Pos.Y*10), 2, c)
        }
    }

    f, err := os.OpenFile(filename + ".png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}


func drawEdgeList(v *Voronoi, edges []EdgeIndex, filename string) {
    var w, h int = 1000, 1000

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // Edges
    // yellow
    //c := color.RGBA{255,255,0,255}
    // blue
    c := color.RGBA{0,100,100,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(2)
    gc.SetStrokeColor(c)

    for i,eI := range edges {

        if eI == EmptyEdge {
            continue
        }

        e := v.edges[eI]

        var tmp HEEdge
        e2   := v.edges[e.ETwin]
        if e != tmp && e2 != tmp {

            edge := Edge{}
            v1   := e.VOrigin
            v2   := e2.VOrigin
            switch {
                // Best case. We have both endpoints. And none of them Infinity ones.
                case v1.Valid() && v2.Valid() && v.vertices[v1].Pos != InfinitePoint && v.vertices[v2].Pos != InfinitePoint:
                    //fmt.Printf("Edge in question: \n\t%v\n\t%v\n\t%v, %v\n", e, e2, v1, v2)
                    edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                // We have the "left" endpoint.
                case v1.Valid() && v.vertices[v1].Pos != InfinitePoint && !v2.Valid():
                    edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                // We have the "right" endpoint.
                case !v1.Valid() && v2.Valid() && v.vertices[v2].Pos != InfinitePoint:
                    edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                // We don't have any endpoints.
                default:
                    //i = i
                    // amplified line. Exceeding all boundaries. Infinite line.
                    edge = createLine(v, EdgeIndex(i), false)
            }



            if edge == (Edge{}) {
                //continue
            }

            //fmt.Printf("From %v|%v ----> %v|%v\n", edge.Pos.X*10., float64(h) - edge.Pos.Y*10., Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)

            gc.MoveTo(edge.Pos.X*10., float64(h) - edge.Pos.Y*10.)
            gc.LineTo(Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)
            gc.FillStroke()
            gc.Close()
        }
    }

    f, err := os.OpenFile(filename + ".png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

func (v *Voronoi)drawFaces() {

    for i,f := range v.faces {
        var edgeIndexList []EdgeIndex

        startEdge := f.EEdge

        edgeIndexList = append(edgeIndexList, startEdge)

        edgeCount := 0
        e := v.edges[startEdge].ENext

        for e != EmptyEdge && e != startEdge {
            edgeIndexList = append(edgeIndexList, e)
            if edgeCount > 50 {
                fmt.Printf("error because of too large edge count?")
                break
            }

            e = v.edges[e].ENext
            edgeCount += 1
        }

        drawEdgeList(v, edgeIndexList, "edges_face_" + fmt.Sprintf("%v", i))

    }

}

func drawDividingChain(chain []ChainElem) {
    var w, h int = 1000, 1000

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // blue
    c := color.RGBA{0,0,255,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(2)

    gc.SetStrokeColor(c)

    empty := Vector{}
    for i,ch := range chain {
        if ch.intersection != empty {
            drawCircle(m, int(ch.intersection.X*10), h-int(ch.intersection.Y*10), 3, c)

            if i > 0 {
                gc.MoveTo(chain[i-1].intersection.X*10., float64(h) - chain[i-1].intersection.Y*10.)
                gc.LineTo(ch.intersection.X*10., float64(h) - ch.intersection.Y*10.)
                gc.FillStroke()
                gc.Close()
            }
        }

    }

    f, err := os.OpenFile("dividing_chain_rec_8.png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

// Verifies that the Voronoi is not corrupted or has miscalculated edges/faces
// or any other unvalid stuff.
func (v *Voronoi) Verify() error {

    emptyF := HEFace{}
    emptyE := HEEdge{}
    emptyV := HEVertex{}

    for i,e := range v.edges {
        if e != emptyE {

            // Every valid edge MUST have a valid Twin edge!
            if e.ETwin == EmptyEdge || v.edges[e.ETwin] == emptyE {
                return errors.New(fmt.Sprintf("Edge %v: %v has an invalid twin edge", i, e))
            }

            // Twins must refer to each other!
            if i != int(v.edges[e.ETwin].ETwin) {
                return errors.New(fmt.Sprintf("Edge %v and his Twin %v don't refer to each other", i, e.ETwin))
            }

            // Check if the origin vertex is valid (if there is one)
            if e.VOrigin != EmptyVertex && e.VOrigin != InfiniteVertex && v.vertices[e.VOrigin] == emptyV {
                return errors.New(fmt.Sprintf("Origin vertex of edge %v: %v is invalid", i, e))
            }

            // Check, if the face is valid
            if e.FFace == EmptyFace || v.faces[e.FFace] == emptyF {
                return errors.New(fmt.Sprintf("Face of edge %v: %v is invalid", i, e))
            }

            // if VOrigin is Empty, it MUST be the referenced edge for the face of the edge!
            if (e.VOrigin == EmptyVertex || e.VOrigin == InfiniteVertex) && i != int(v.faces[e.FFace].EEdge) {
                return errors.New(fmt.Sprintf("A first edge %v: %v must be referenced as first edge by the corresponding face %v!", e, i, e.FFace))
            }

            // Check, if the next edge is valid
            if e.ENext != EmptyEdge && v.edges[e.ENext] == emptyE {
                return errors.New(fmt.Sprintf("Next edge for edge %v: %v is invalid", i, e))
            }

        }
    }

    for i,f := range v.faces {
        if f != emptyF {

            // Check for valid reference point
            if f.ReferencePoint == InfinitePoint {
                return errors.New(fmt.Sprintf("Face %v: %v has infinite reference point", i, f))
            }

            // Check for valid EEdge
            if f.EEdge == EmptyEdge || v.edges[f.EEdge] == emptyE {
                return errors.New(fmt.Sprintf("Face %v: %v points to invalid edge", i, f))
            }

            // If the edge is the first one (no/infinite start-vertex), it's all good.
            // Otherwise check, that the edges actually go all the way around the face!
            if v.edges[f.EEdge].VOrigin != EmptyVertex && v.edges[f.EEdge].VOrigin != InfiniteVertex {
                startEdge := f.EEdge
                edgeCount := 0
                e := v.edges[startEdge].ENext

                for e != EmptyEdge && e != startEdge {

                    if int(v.edges[e].FFace) != i {
                        return errors.New(fmt.Sprintf("Edge %v of the face %v: %v does not point to the right face!", e, i, f))
                    }

                    if edgeCount > 50 {
                        return errors.New(fmt.Sprintf("Looping around the edges of the face %v: %v most likely does not stop (infinite loop)", i, f))
                    }

                    e = v.edges[e].ENext
                    edgeCount += 1
                }
            }
        }
    }

    return nil
}

////////////////////////////////////////////////////////////////////////
//  Calculate a Voronoi Tessellation
////////////////////////////////////////////////////////////////////////

// Calculates a convex hull for the given voronoi. Runs in O(n).
// The convex hull consists of a list of Reference points.
// As there is an exactly 1:1 representation of Reference points and HEFaces,
// a list of face indices is returned in Order (counter-clockwise) of the convex hull.
func (v *Voronoi)ConvexHull(fEntry VoronoiEntryFace) ConvexHull {

    startFace := FaceIndex(fEntry)
    convexHullList := []FaceIndex{FaceIndex(startFace)}

    // If the voronoi consists of just one refPoint, there is no edge.
    if v.faces[startFace].EEdge != EmptyEdge {

        nextFace := v.edges[v.edges[v.faces[startFace].EEdge].ETwin].FFace

        convexHullList = append(convexHullList, nextFace)

        for {
            nextFace = v.edges[v.edges[v.faces[nextFace].EEdge].ETwin].FFace

            // Once all around.
            if nextFace == startFace {
                return convexHullList
            }

            lastFace := convexHullList[len(convexHullList)-2]
            // Special case of linear voronoi.
            if nextFace == lastFace {
                return convexHullList
            }

            convexHullList = append(convexHullList, nextFace)
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
// Returns the index of the convex hull list called upon.
func (ch ConvexHull) bestFaceIndex(v *Voronoi, isBetter func(v1, v2 Vector) bool) int {
    bestFace := 0
    //fmt.Println ("  bestFace:", ch)
    for i,_ := range ch {
        if isBetter(v.faces[ch[i]].ReferencePoint, v.faces[ch[bestFace]].ReferencePoint) {
            bestFace = i
        }
    }
    return bestFace
}

// Creates a line from the HEEdge. Depeding on its state from the edge or the TmpEdge.
func createLine (v *Voronoi, e EdgeIndex, amplified bool) Edge {



    switch {
        case !v.edges[e].VOrigin.Valid() && v.edges[v.edges[e].ETwin].VOrigin.Valid():
            edge := v.edges[e].TmpEdge
            edge.Dir = Mult(edge.Dir, 50)
            return Amplify(edge, 50)
        case v.edges[e].VOrigin.Valid() && !v.edges[v.edges[e].ETwin].VOrigin.Valid():
            edge := v.edges[e].TmpEdge
            edge.Dir = Mult(edge.Dir, 50)
            return Amplify(edge, 50)
        case !v.edges[e].VOrigin.Valid() && !v.edges[v.edges[e].ETwin].VOrigin.Valid():
             return Amplify(v.edges[e].TmpEdge, 50.0)
        default:
            //fmt.Println("No problems here!!!!!!!!!!!!!!!!")
            return Edge {
                Pos: v.vertices[v.edges[e].VOrigin].Pos,
                Dir: Sub(v.vertices[v.edges[v.edges[e].ETwin].VOrigin].Pos, v.vertices[v.edges[e].VOrigin].Pos),
            }
    }

    /*
    if !v.edges[e].VOrigin.Valid() || !v.edges[v.edges[e].ETwin].VOrigin.Valid() {
        if amplified {
            return Amplify(v.edges[e].TmpEdge, 50.0)
        } else {
            return v.edges[e].TmpEdge
        }
    } else {
        //fmt.Println("No problems here!!!!!!!!!!!!!!!!")
        return Edge {
            Pos: v.vertices[v.edges[e].VOrigin].Pos,
            Dir: Sub(v.vertices[v.edges[v.edges[e].ETwin].VOrigin].Pos, v.vertices[v.edges[e].VOrigin].Pos),
        }
    }
    */
}

// Calculates the highest intersection of the given bisector and any edge of the
// given face. With the restriction, that it can't be 'lastEdge'.
func calcHighestIntersection(v *Voronoi, bisector Edge, face FaceIndex, lastEdge EdgeIndex, lastVertex Vector) (EdgeIndex, Vector) {

    //fmt.Println("begin--------------------------------------------------")

    bestEdge := EmptyEdge
    bestIntersection := Vector{}

    minY := -100000000000.0
    maxY := -minY
    if lastVertex != (Vector{}) {
        maxY = lastVertex.Y
    }

    edge := v.faces[face].EEdge
    veryFirstEdge := edge

    first := true


    // Find be highest/best intersection of the bisector and edge.
    for edge != EmptyEdge && (first || edge != veryFirstEdge) {

        intersects, location := LineIntersection4(bisector, createLine(v, edge, true))
        //fmt.Printf("Found an intersection: %v, %v\n", intersects, location)

        if (g_recursions == 8) {
            g_edgeIndexList = append(g_edgeIndexList, edge)
        }

        if intersects &&
           edge != lastEdge &&
           location.Y < maxY &&
           location.Y > minY &&
           !Equal(location, lastVertex) &&
           !Equal(location, bestIntersection) {
            //fmt.Printf("--> Better than the last! %v\n", location)
            bestEdge = edge
            bestIntersection = location
            minY = bestIntersection.Y


        }

        first = false
        edge = v.edges[edge].ENext
    }
    //fmt.Println("end----------------------------------------------------")
    return bestEdge, bestIntersection
}

// Calculates the upper or lower common support line for the two given convex hulls.
// Depending on the betterSide and betterPoint-function, it could potentially also calculate
// left/right common support lines, depending on the divide-and-conquer approach of the main algorithm.
func commonSupportLine(v *Voronoi, h1, h2 ConvexHull, betterSideL func(v1, v2 Vector) bool, betterSideR func(v1, v2 Vector) bool, betterPoint func(v1, v2, test Vector) bool) (FaceIndex, FaceIndex) {

    ai := h1.bestFaceIndex(v, betterSideL)
    bi := h2.bestFaceIndex(v, betterSideR)

    finished := false

    for !finished {
        finished = true

        nextAi := (ai-1)%len(h1)
        if nextAi < 0 {
            nextAi = len(h1)-1
        }

        // iterating through the convex hull points...
        if len(h1) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[nextAi]].ReferencePoint) {
            ai = nextAi
            finished = false
        }
        // iterating through the convex hull points...
        if len(h2) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[(bi+1)%len(h2)]].ReferencePoint) {
            bi = (bi+1)%len(h2)
            finished = false
        }

        // We try both directions every time! Because of some special cases, where the first point for h1 is the correct one, but for the very
        // first step (maybe later?), the next point is better for the moment. It then never goes back to the actual globally better point.
        // By checking both directions, we can actually get back. And there should be no downsides, as there is only ever one direction to go!

        // iterating through the convex hull points...
        if len(h1) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[(ai+1)%len(h1)]].ReferencePoint) {
            ai = (ai+1)%len(h1)
            finished = false
        }

        nextBi := (bi-1)%len(h2)
        if nextBi < 0 {
            nextBi = len(h2)-1
        }

        // iterating through the convex hull points...
        if len(h2) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[nextBi]].ReferencePoint) {
            bi = nextBi
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
        return test.X >= v1.X && test.X <= v2.X
    }
    return false
}

//
// Extract the dividing chain between the two voronois without actually
// merging them! This will be done in an additional step.
//
func (v *Voronoi)extractDividingChain(left, right VoronoiEntryFace) []ChainElem {
    var dividingChain []ChainElem

    fmt.Printf("dividing chain for recursion: %v\n", g_recursions)

    h1 := v.ConvexHull(left)
    h2 := v.ConvexHull(right)

    isRight := func(v1, v2 Vector) bool {
        return v1.X > v2.X
    }
    isLeft := func(v1, v2 Vector) bool {
        return v1.X < v2.X
    }
    // Upper common support line!
    p, q := commonSupportLine(v, h1, h2, isLeft, isRight, isBetterUp)
    // Lower common support line!
    h1Down, h2Down := commonSupportLine(v,  h1, h2, isRight, isLeft, isBetterDown)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := Vector{}

    bisector := PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
    bisector = Amplify(bisector, 100.0)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge, lastVertex)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge, lastVertex)

        switch {

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case edgeP == EmptyEdge && edgeQ == EmptyEdge:
                dividingChain = append(dividingChain, ChainElem{Vector{}, EmptyEdge, EmptyEdge, p, q, bisector})

            // Equal Intersection with both edges. So 4 edges meet.
            case edgeQ != EmptyEdge && edgeP != EmptyEdge && Equal(locationP, locationQ):
                // This could result in kind of an invalid Delaunay triangulations. At least regarding a triangle.
                // This should now work for situations, where 4 edges meet. I have to re-examine for even more edges meeting...

                dividingChain = append(dividingChain, ChainElem{locationP, edgeP, edgeQ, p, q, bisector})

                lastVertex   = locationP
                lastPEdge    = edgeP
                lastQEdge    = edgeQ

                p = v.edges[v.edges[edgeP].ETwin].FFace
                q = v.edges[v.edges[edgeQ].ETwin].FFace

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && (edgeQ == EmptyEdge || locationP.Y >= locationQ.Y):

                dividingChain = append(dividingChain, ChainElem{locationP, edgeP, EmptyEdge, p, q, bisector})

                lastVertex   = locationP
                lastPEdge    = edgeP
                lastQEdge    = EmptyEdge

                p = v.edges[v.edges[edgeP].ETwin].FFace

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && (edgeP == EmptyEdge || locationQ.Y >= locationP.Y):

                dividingChain = append(dividingChain, ChainElem{locationQ, EmptyEdge, edgeQ, p, q, bisector})

                lastVertex   = locationQ
                lastQEdge    = edgeQ
                lastPEdge    = EmptyEdge

                q = v.edges[v.edges[edgeQ].ETwin].FFace


            default:
                fmt.Println("default")
                lastMerge = true
                v.pprint()
                os.Exit(0)
        }

        if lastMerge {
            break
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
        bisector = Amplify(bisector, 100.0)

    }

    if g_recursions == 7 {
        //drawDividingChain(dividingChain)
    }

    return dividingChain
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func (v *Voronoi)mergeVoronoi(left, right VoronoiEntryFace) VoronoiEntryFace {
    var emptyV HEVertex
    var emptyE HEEdge

    g_recursions += 1;

    nextPEdge    := EmptyEdge
    lastDownEdge := EmptyEdge
    lastVertex   := EmptyVertex

    lastQRemove  := EmptyEdge
    lastQVertex  := EmptyVertex

    potFirstEdgeOfQ := EmptyEdge
    lastQFace    := EmptyFace

    qFirstContact := false
    //specialQEdge := EmptyEdge
    var specialQEdges []EdgeIndex

    chainList := v.extractDividingChain(left, right)

    // Iterating through the dividing chain and actually merge the voronois.
    for i,chain := range chainList {
        fmt.Printf("i == %v\n", i)
        heVertex := EmptyVertex
        if chain.intersection != (Vector{}) {

            if chain.intersection == InfinitePoint {
                fmt.Println("INFINITY")
                heVertex = InfiniteVertex
            } else {
                heVertex = v.createVertex(chain.intersection)
            }
        }

        otherWayBisector := chain.bisector
        otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

        heEdgeUp   := v.createEdge(heVertex,   EmptyEdge, nextPEdge,   chain.p, otherWayBisector)
        heEdgeDown := v.createEdge(lastVertex, heEdgeUp,  chain.edgeQ, chain.q, chain.bisector)
        v.edges[heEdgeUp].ETwin = heEdgeDown

        if lastDownEdge != EmptyEdge {
            v.edges[lastDownEdge].ENext = heEdgeDown
        }

        // For the very first dividing chain element holds, that the face q has to update
        // its edge reference to this down-edge. No matter what.
        if i == 0 {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if i == len(chainList)-1 {
            v.faces[chain.p].EEdge = heEdgeUp
        }

        // If edgeQ is the current edge for the face q and we have an intersection with
        // edgeQ, the new first edge for the face q must be edgeQ!
        // Simple case!
        if chain.edgeQ != EmptyEdge && v.faces[chain.q].EEdge == chain.edgeQ && chain.intersection != InfinitePoint {
            fmt.Printf("are we here? %v, %v, %v\n", chain.q, heEdgeDown, chain.edgeQ)
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if v.faces[chain.q].EEdge == EmptyEdge || (lastDownEdge == EmptyEdge && heVertex != InfiniteVertex) {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if v.faces[chain.p].EEdge == EmptyEdge || !heVertex.Valid() {
            v.faces[chain.p].EEdge = heEdgeUp
        }

        // If we just enter a q-face polygon for the first time!
        // Then we save the heEdgeDown because it is a (very) good candidate for the first edge of q in times of crisis.
        if i > 0 && chainList[i-1].edgeQ == EmptyEdge && chain.edgeQ != EmptyEdge {
            potFirstEdgeOfQ = heEdgeDown
            lastQFace = chain.q
            qFirstContact = true
            fmt.Printf("q == %v, potFirstEdgeOfQ == %v\n", lastQFace, potFirstEdgeOfQ)
        }

        // Removing stuff!
        if chain.edgeP != EmptyEdge {
            // Delete vertices that overlap with the other side (to the left)!
            vertex := v.edges[v.edges[chain.edgeP].ETwin].VOrigin

            if vertex.Valid() && v.vertices[vertex] != emptyV {
                fmt.Printf("    ----> Found a vertex that has do be deleted (P) (v-index: %v - %v)\n", vertex, v.vertices[vertex].Pos)

                if v.edges[v.edges[chain.edgeP].ENext] != emptyE {
                    if v.faces[v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin].FFace].EEdge == v.edges[v.edges[chain.edgeP].ENext].ETwin {
                        // The edge is referenced as first edge. Might need some extra work here...
                        v.faces[v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin].FFace].EEdge = v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin].ENext
                    }

                    v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin] = emptyE

                }

                if v.edges[v.edges[chain.edgeP].ENext].ENext != EmptyEdge {
                    v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin].VOrigin = EmptyVertex
                    v.vertices[v.edges[v.edges[v.edges[chain.edgeP].ENext].ENext].VOrigin] = emptyV
                    v.edges[v.edges[v.edges[v.edges[chain.edgeP].ENext].ENext].ETwin] = emptyE
                    v.edges[v.edges[v.edges[chain.edgeP].ENext].ENext] = emptyE
                    v.edges[v.edges[chain.edgeP].ENext].ENext = EmptyEdge
                }

                v.edges[v.edges[v.edges[chain.edgeP].ENext].ETwin] = emptyE

                // Delete the following edge of edgeP
                v.edges[v.edges[chain.edgeP].ENext] = emptyE
                v.edges[chain.edgeP].ENext = EmptyEdge

                // Delete the vertex.
                v.vertices[v.edges[v.edges[chain.edgeP].ETwin].VOrigin] = emptyV
                v.edges[v.edges[chain.edgeP].ETwin].VOrigin = EmptyVertex
            }
        }

        // Removing Stuff
        if chain.edgeQ != EmptyEdge {
            // Delete vertices that overlap with the other side (to the right)!
            vertex := v.edges[chain.edgeQ].VOrigin
            if vertex.Valid() && v.vertices[vertex] != emptyV {
                fmt.Printf("    ----> Found a vertex that has do be deleted (Q:%v) (v-index: %v)\n", chain.q, vertex)
                fmt.Printf("    q.EEdge: %v\n", v.faces[chain.q].EEdge)
                test := true

                if test {

                    // reset old states
                    if lastQRemove == chain.edgeQ {
                        fmt.Printf("lastQRemove set to EmptyEdge\n")
                        lastQRemove = EmptyEdge
                    }

                    // THIS
                    //   \  /x
                    //    \/
                    //    | /x
                    //     /
                    //     \
                    //      \x
                    // Save a reference to THIS special edge. It will be removed in the last step.
                    if qFirstContact {
                        //specialQEdge = v.edges[v.edges[v.edges[v.edges[chain.edgeQ].ETwin].ENext].ETwin].ENext
                        specialQEdges = append(specialQEdges, v.edges[v.edges[v.edges[v.edges[chain.edgeQ].ETwin].ENext].ETwin].ENext)
                    }

                    // remove edges
                    if lastQRemove != EmptyEdge {
                        fmt.Printf("Condition: %v || %v\n", v.edges[lastQRemove].ETwin == v.faces[chain.q].EEdge, lastQRemove == v.faces[chain.q].EEdge)
                        // well, here it is, our crisis.
                        if v.edges[lastQRemove].ETwin == v.faces[chain.q].EEdge || lastQRemove == v.faces[chain.q].EEdge {
                            fmt.Printf("Nearly lost the reference to the first face of q.\n")
                            v.faces[chain.q].EEdge = potFirstEdgeOfQ
                        }
                        fmt.Printf("DELETE E: %v\n", v.edges[lastQRemove].ETwin)
                        v.edges[v.edges[lastQRemove].ETwin] = emptyE
                        fmt.Printf("DELETE E: %v\n", lastQRemove)
                        v.edges[lastQRemove] = emptyE
                    }
                    // remove old vertex anywhay
                    if lastQVertex != EmptyVertex {
                        fmt.Printf("DELETE V: %v: %v\n", lastQVertex, v.vertices[lastQVertex])
                        v.vertices[lastQVertex] = emptyV
                    }

                    // set new states
                    lastQRemove = v.edges[v.edges[chain.edgeQ].ETwin].ENext

                    if vertex != lastQVertex {
                        lastQVertex = vertex
                    } else {
                        lastQVertex = EmptyVertex
                    }

                } else {

                    edgeRef := v.edges[v.edges[v.edges[chain.edgeQ].ETwin].ENext].ETwin
                    edgeToBeDeleted := v.edges[v.edges[edgeRef].ENext].ETwin

                    fmt.Printf("SETFACE: %v --> %v\n", chain.q, heEdgeDown)
                    v.faces[chain.q].EEdge = heEdgeDown

                    // Delete Edge
                    //fmt.Printf("q edgeRef: %v\n", edgeRef)
                    //v.edges[edgeRef].ENext = EmptyEdge

                    fmt.Printf("DELETE E: %v\n", v.edges[edgeToBeDeleted].ETwin)
                    v.edges[v.edges[edgeToBeDeleted].ETwin] = emptyE
                    // Delete Edge
                    fmt.Printf("DELETE E: %v\n", edgeToBeDeleted)
                    v.edges[edgeToBeDeleted] = emptyE
                    // Delete Vertex
                    fmt.Printf("DELETE V: %v: %v\n", vertex, v.vertices[vertex])
                    v.vertices[vertex] = emptyV

                }



            }
        }

        // Setting States
        if chain.edgeP != EmptyEdge {
            v.edges[chain.edgeP].ENext = heEdgeUp
            v.edges[v.edges[chain.edgeP].ETwin].VOrigin = heVertex
            nextPEdge = v.edges[chain.edgeP].ETwin
        } else {
            nextPEdge = heEdgeUp
        }

        // Setting States
        if chain.edgeQ != EmptyEdge {
            v.edges[chain.edgeQ].VOrigin = heVertex
            lastDownEdge = v.edges[chain.edgeQ].ETwin
        } else {
            lastDownEdge = heEdgeDown
        }

        lastVertex = heVertex
        qFirstContact = false

    }

    // remove edges
    if lastQRemove != EmptyEdge && v.edges[lastQRemove] != emptyE {
        // well, here it is, our crisis.
        //lastQ := v.edges[v.edges[lastQRemove].ETwin].FFace
        lastQ := lastQFace

        //fmt.Printf("is it true? %v\n", lastQ == lastQFace)

        if v.edges[lastQRemove].ETwin == v.faces[lastQ].EEdge || lastQRemove == v.faces[lastQ].EEdge {
            //fmt.Printf("Nearly lost the reference to the first face of q.\n")
            fmt.Printf("LAST SETFACE: %v --> %v\n", lastQ, potFirstEdgeOfQ)
            v.faces[lastQ].EEdge = potFirstEdgeOfQ
        }

        endangeredFace := v.edges[v.edges[lastQRemove].ETwin].FFace
        fmt.Printf("face: %v\n", endangeredFace)

        fmt.Printf("LAST DELETE E: %v\n", v.edges[lastQRemove].ETwin)
        v.edges[v.edges[lastQRemove].ETwin] = emptyE
        fmt.Printf("LAST DELETE E: %v\n", lastQRemove)
        v.edges[lastQRemove] = emptyE

        if v.edges[v.faces[endangeredFace].EEdge] == emptyE {
            v.faces[endangeredFace].EEdge = lastDownEdge
        }
    }
    // remove old vertex anywhay
    if lastQVertex != EmptyVertex {
        fmt.Printf("LAST DELETE V: %v\n", lastQVertex)
        v.vertices[lastQVertex] = emptyV
    }
    // potentially remove the special edge of q
    for _,e := range specialQEdges  {
        if e != EmptyEdge && v.edges[e] != emptyE {
            fmt.Printf("SPECIAL DELETE E: %v\n", v.edges[e].ETwin)
            v.edges[v.edges[e].ETwin] = emptyE
            fmt.Printf("SPECIAL DELETE E: %v\n", e)
            v.edges[e] = emptyE
        }
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
        vertices:           make([]HEVertex, 2*n-5 + 3 + 100),
        firstFreeVertexPos: 0,
        edges:              make([]HEEdge, 2*(3*n-6) + 6 + 100),
        firstFreeEdgePos:   0,
        faces:              make([]HEFace, n),
        firstFreeFacePos:   0,
    }

    v.divideAndConquer(pointList)

    e := v.Verify()

    if e != nil {
        fmt.Printf("VERIFICATION_ERROR: '%v'\n", e)
    }

    return v
}

//
// Normal test distribution of points. No special cases for Voronoi generation.
//
func testNormal01() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 50., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_01", true)

}

func testNormal02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 30., 0})
    pointList = append(pointList, Vector{80., 20., 0})
    pointList = append(pointList, Vector{55., 40., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_02", true)
}

func testNormal03() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_03\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{45., 45., 0})
    pointList = append(pointList, Vector{50., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_03", true)

}

//
// Test cases where 4 lines meet and p/q have the same intersection point.
//
func testEqualIntersection01() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_intersection_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{28, 30, 0})
    pointList = append(pointList, Vector{30, 20, 0})
    pointList = append(pointList, Vector{30, 40, 0})

    pointList = append(pointList, Vector{40, 20, 0})
    pointList = append(pointList, Vector{40, 40, 0})
    pointList = append(pointList, Vector{42, 30, 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_equal_intersection_01", true)
}

func testEqualIntersection02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_intersection_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 50., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{60., 20., 0})
    pointList = append(pointList, Vector{70., 50., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_equal_intersection_02", true)
}

//
// Test cases with linear dependent points or lines
//
func testLinearDepentence01() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20, 10, 0})
    pointList = append(pointList, Vector{30, 10, 0})
    pointList = append(pointList, Vector{40, 10, 0})
    pointList = append(pointList, Vector{50, 10, 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_01", true)
}

func testLinearDepentence02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_02", true)
}

func testLinearDepentence03() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_03\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20., 10., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 10., 0})
    pointList = append(pointList, Vector{80., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_03", true)
}

func testLinearDepentence04() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_04\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_04", true)
}

//
// Test cases that cause infinite loops or unknown problems.
//

func testUnknownProblem01() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_01\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 15., 0})
    pointList = append(pointList, Vector{15., 65., 0})
    pointList = append(pointList, Vector{20., 40., 0})

    pointList = append(pointList, Vector{40.,  5., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 60., 0})
    // Adding this point results in a wrong voronoi! Really messed up.
    pointList = append(pointList, Vector{90., 10., 0})

    for i,_ := range pointList {
        pointList[i] = Mult(pointList[i], 0.5)
        pointList[i] = Add(pointList[i], Vector{25,25,0})
    }

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_01", true)
}

func testUnknownProblem02() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_02\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})

    pointList = append(pointList, Vector{25., 40., 0})

    pointList = append(pointList, Vector{30., 10., 0})

    pointList = append(pointList, Vector{40., 30., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_02", true)
}

func testUnknownProblem03() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_03\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 15., 0})
    pointList = append(pointList, Vector{15., 65., 0})
    pointList = append(pointList, Vector{20., 40., 0})

    pointList = append(pointList, Vector{40.,  5., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 60., 0})

    // Infinite loop by adding this point.
    pointList = append(pointList, Vector{70., 70., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_03", true)
}

func testUnknownProblem04() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_04\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{26.140400405282577, 61.621512929924580, 0})
    pointList = append(pointList, Vector{30.777003166789225, 25.899889187518045, 0})
    pointList = append(pointList, Vector{38.049244986921730, 50.778804685015940, 0})
    pointList = append(pointList, Vector{39.182654173412490, 53.260689369108890, 0})
    pointList = append(pointList, Vector{56.653551080486500, 67.211346445470530, 0})

    pointList = append(pointList, Vector{61.682906469053240, 71.811352995110340, 0})
    pointList = append(pointList, Vector{61.787090244227150, 57.606346468267920, 0})
    pointList = append(pointList, Vector{69.598081344490680, 56.223764727531960, 0})
    pointList = append(pointList, Vector{71.996875920571910, 27.832275576849792, 0})
    pointList = append(pointList, Vector{72.329734102573130, 63.247795405728404, 0})


    if false {
        for i,_ := range pointList {
            //pointList[i] = Mult(pointList[i], 0.5)
            pointList[i] = Add(pointList[i], Vector{-20,0,0})
        }
    }

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_04", true)
}

func testUnknownProblem05() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_05\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409032579303917
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_05", true)
}

func testUnknownProblem06() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_06\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409509781449406
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }
    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_06", true)
}

func testUnknownProblem07() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_07\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409694089411269
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }
    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_07", true)
}

func testUnknownProblem08() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_08\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409796721354066
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_08", true)
}

func testUnknownProblem09() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_09\n")
    fmt.Printf("===========================\n")
    count := 5
    var seed int64 = 1482476817056249360
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        //v := Vector{r.Float64()*150., r.Float64()*150., 0}
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[2:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_09", false)
}

func testUnknownProblem10() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_10\n")
    fmt.Printf("===========================\n")
    count := 5
    var seed int64 = 1482477466512747360
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_10", true)
}

func testUnknownProblem11() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_11\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409032579303917
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)
    pointList = pointList[:8]

    //pointList = pointList[:4]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_11", true)
}

func testUnknownProblem12() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_12\n")
    fmt.Printf("===========================\n")
    count := 9
    var seed int64 = 1483370038194119290
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)
    pointList = pointList[1:]

    //pointList = pointList[4:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_12", true)
}

func testUnknownProblemSeed(seed int64, count int) {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_seed_%v\n", strconv.FormatInt(seed, 10))
    fmt.Printf("===========================\n")
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    pointList = pointList[len(pointList)/2:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_seed_" + strconv.FormatInt(seed, 10), true)

    //v.drawFaces()

}

//
// Test cases with random points
//

func testRandom(count int) {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_random_%v\n", count)
    fmt.Printf("===========================\n")
    seed := time.Now().UTC().UnixNano()
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    fmt.Printf("Seed: %v\n", strconv.FormatInt(seed, 10))
    fmt.Printf("Random poins:\n\t")
    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
        fmt.Printf("%v, ", v)
    }
    fmt.Printf("\n\n")

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_random_" + strconv.FormatInt(seed, 10), true)
}


func main() {

    working := false

    if working {
        testNormal01()
        testNormal02()
        testNormal03()
        // crash
        //testEqualIntersection01()
        testEqualIntersection02()
        testLinearDepentence04()
        testUnknownProblem03()
        testUnknownProblem09()
        testUnknownProblem10()
        testUnknownProblemSeed(1489941049442429888, 5)
        testUnknownProblem01()
        testUnknownProblem08()
        testUnknownProblem11()
        testUnknownProblem06()
        // wrong
        //testUnknownProblem07()
        // wrong
        //testUnknownProblem12()
        testUnknownProblem04()
        // wrong
        //testUnknownProblemSeed(1483369884537650258, 20)
        testUnknownProblemSeed(1483370089898481236, 15)
        testUnknownProblem05()
        testUnknownProblemSeed(1483370089898481236, 15)
        testUnknownProblem02()
    }

    toBeVerified := false

    if toBeVerified {

        testLinearDepentence01()
        // works.
        // A first edge {-1 0 6 1 {{35 510 0} {0 -1000 0}}}: 1 must be referenced as first edge by the corresponding face 1 !
        fmt.Println("Test: testLinearDepentence02")
        testLinearDepentence02()

        // works.
        // A first edge {-2 1 -1 1 {{35 510 0} {-0 1000 -0}}}: 0 must be referenced as first edge by the corresponding face 1 !
        fmt.Println("Test: testLinearDepentence03")
        testLinearDepentence03()

        // works.
        for i := 0; i < 20; i++ {
            fmt.Println("Test: testRandom_", i)
            testRandom(5)
        }

        // Edge 1: {6 0 26 2 {{241.90357878309425 103.13635491102873 0} {-427.54352793789894 -125.38564346969068 0}}} has an invalid twin edge
        testUnknownProblemSeed(1483370150842201370, 15)

        // crashes...
        testUnknownProblemSeed(1483370130545841965, 15)

    }

    crashes := false

    if crashes {
        // Both crashes can be resolved by increasing the static size for vertices
        // and edges.
        // Idea: Add an additional list for vertices, edges and faces (?) that contain a list
        // of free additional positions for new vertices/edges.
        // On those lists, only prepend-operations will be done [O(1)].
        // If len(list) > 0, we just take the first additional free position. That way, we
        // never get over the max list size and do not get additional time other than O(1).

    }

    infLoop := false

    if infLoop {

    }

    test := true

    if test {
        testUnknownProblemSeed(1483370150842201370, 15)
    }

}




