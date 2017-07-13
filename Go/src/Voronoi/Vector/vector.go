package vector

import (
    //"fmt"
    "math"
    //"reflect"
)

const (
    EPS = 0.00001
)

var InfinitePoint = Vector{-100000000,-100000000, -100000000}

// Because of laziness, a Vector could also just be a point. Depending on context.
type Vector struct {

    // Potentially changing to float128 --> https://gist.github.com/grd/4050062
    // Or math/big --> https://golang.org/pkg/math/big/

    X,Y,Z   float64
}

type Edge struct {
    Pos             Vector
    Dir             Vector
}

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

// Just makes the edge a lot bigger according to s.
func (e *Edge) Amplify(s float64) {

    //e.Dir = Mult(e.Dir,s)
    //e.Pos = Sub(e.Pos, Mult(e.Dir, 0.5))

    e.Dir.Mult(s)
    e.Pos.Sub(e.Dir)
    e.Dir.Mult(2)
}

func Amplify(e Edge, s float64) Edge {
    return Edge{
        Dir: Mult(e.Dir,s),
        Pos: Add(e.Pos, Mult(e.Dir, -s/2.0)),
    }
}

func (e Edge) Copy() Edge {
    return Edge{Pos: e.Pos.Copy(), Dir: e.Dir.Copy()}
}

func Equal(v1, v2 Vector) bool {
    return  math.Abs(v1.X-v2.X) <= EPS &&
            math.Abs(v1.Y-v2.Y) <= EPS &&
            math.Abs(v1.Z-v2.Z) <= EPS
}

func Add(v1, v2 Vector) Vector {
    return Vector{X: v1.X+v2.X, Y: v1.Y+v2.Y, Z: v1.Z+v2.Z}
}
func Mult(v1 Vector, s float64) Vector {
    return Vector{X: v1.X*s, Y: v1.Y*s, Z: v1.Z*s}
}
func Sub(v1, v2 Vector) Vector {
    return Vector{X: v2.X-v1.X, Y: v2.Y-v1.Y, Z: v2.Z-v1.Z}
}

func (v Vector) Copy() Vector {
    return Vector{X:v.X, Y:v.Y, Z:v.Z}
}

func (v *Vector) Add(v1 Vector) {
    v.X += v1.X
    v.Y += v1.Y
    v.Z += v1.Z
}

func (v *Vector) Sub(v1 Vector) {
    v.X -= v1.X
    v.Y -= v1.Y
    v.Z -= v1.Z
}

func (v *Vector) Mult(s float64) {
    v.X *= s
    v.Y *= s
    v.Z *= s
}

// Calculates v2 - v1.
func Diff (v1, v2 Vector) Vector {
    return Vector {
        X: v2.X - v1.X,
        Y: v2.Y - v1.Y,
        Z: v2.Z - v1.Z,
    }
}

// Calculates the cross product: v1 x v2
func Cross(v1, v2 Vector) Vector {
    return Vector {
        X: v1.Y*v2.Z - v1.Z*v2.Y,
        Y: v1.Z*v2.X - v1.X*v2.Z,
        Z: v1.X*v2.Y - v1.Y*v2.X,
    }
}

// Calculates the middle point between p1 and p2.
func MiddlePoint(p1, p2 Vector) Vector {
    d := Diff(p1, p2)
    return Vector {
        X: p1.X + d.X/2.0,
        Y: p1.Y + d.Y/2.0,
        Z: p1.Z + d.Z/2.0,
    }
}

// Calculates an edge representation of one perpendicular bisector between p1 and p2
func PerpendicularBisector(p1, p2 Vector) Edge {
    return Edge {
        Pos: MiddlePoint(p1, p2),
        Dir: Cross(Diff(p1, p2), Vector{0,0,1}),
    }
}

func SideOfLine(v1, v2, test Vector) float64 {
    return (v2.X - v1.X) * (test.Y - v1.Y) - (v2.Y - v1.Y) * (test.X - v1.X)
}

// To determine the higher/lower common support lines of the convex hulls.
func IsLeft2D(v1, v2, test Vector) bool {
    return SideOfLine(v1, v2, test) > 0
}
func IsRight2D(v1, v2, test Vector) bool {
    return SideOfLine(v1, v2, test) < 0
}

func LineIntersection4(e1 Edge, e2 Edge) (bool, Vector) {

    e1E := Add(e1.Pos, e1.Dir)
    //e2E := Add(e2.Pos, e2.Dir)

    det := e1.Dir.Y * e2.Dir.X - e1.Dir.X * e2.Dir.Y
    s1 := (e1E.X-e1.Pos.X) * (e2.Pos.Y-e1.Pos.Y) - (e1E.Y-e1.Pos.Y) * (e2.Pos.X-e1.Pos.X)
    //s2 := (e2E.X-e2.Pos.X) * (e2.Pos.Y-e1.Pos.Y) - (e2E.Y-e2.Pos.Y) * (e2.Pos.X-e1.Pos.X)

    // Lines are colinear.
    if math.Abs(det) <= EPS {
        //return false, Vector{}
        return true, InfinitePoint
    }

    // Intersection is outside of either line, if (s1||s2)/det is outside of 0..1.
    return true, Add(e2.Pos, Mult(e2.Dir, s1/det))
}

// True, if p lies on l.
func PointOnLine(l Edge, p Vector) bool {

    p1 := l.Pos
    p2 := Add(l.Pos, l.Dir)

    dxl := float64(p2.X - p1.X)
    dyl := float64(p2.Y - p1.Y)

    if math.Abs(dxl) >= math.Abs(dyl) {
        if dxl > 0 {
            return p1.X <= p.X && p.X <= p2.X
        } else {
            return p2.X <= p.X && p.X <= p1.X
        }
    } else {
        if dyl > 0 {
            return p1.Y <= p.Y && p.Y <= p2.Y
        } else {
            return p2.Y <= p.Y && p.Y <= p1.Y
        }
    }
}

