package vector

import (
    "fmt"
    "github.com/paulmach/go.geo"
)

// Because of laziness, a Vector could also just be a point. Depending on context.
type Vector struct {
    X,Y,Z   float32
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
func (e *Edge) Amplify(s float32) {
    e.Dir.Mult(s)
    e.Pos.Sub(e.Dir)
    e.Dir.Mult(2)
}

func (e Edge) Copy() Edge {
    return Edge{Pos: e.Pos.Copy(), Dir: e.Dir.Copy()}
}


func Add(v1, v2 Vector) Vector {
    return Vector{X: v1.X+v2.X, Y: v1.Y+v2.Y, Z: v1.Z+v2.Z}
}
func Mult(v1 Vector, s float32) Vector {
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

func (v *Vector) Mult(s float32) {
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

// Determinante of v1 and v2, ignoring the Z-component!
// Honestly, I don't really know, what this does^^
func Det2D(v1, v2 Vector) float32 {
    return v1.X*v2.Y - v1.Y*v2.X
}

// Calculates the intersection point of v1 and v2, if it exists.
// http://stackoverflow.com/questions/20677795/find-the-point-of-intersecting-lines
func LineIntersection(e1 Edge, e2 Edge) Vector {

    xdiff := Vector{X: e1.Dir.X, Y: e2.Dir.X}
    ydiff := Vector{X: e1.Dir.Y, Y: e2.Dir.Y}

    div := Det2D(xdiff, ydiff)
    if div <= 0.0001 {
        fmt.Println("Lines do NOT intersect!")
        return Vector{}
    }

    d := Vector{X: Det2D(e1.Pos, Add(e1.Pos, e1.Dir)), Y: Det2D(e2.Pos, Add(e2.Pos, e2.Dir))}

    return Vector{
        X: Det2D(d, xdiff) / div,
        Y: Det2D(d, ydiff) / div,
        Z: 0.0,
    }
}

/**
inline bool lines_intersect_2d(Vector2 const& p0, Vector2 const& p1, Vector2 const& p2, Vector2 const& p3, Vector2* i const = 0) {
    Vector2 const s1 = p1 - p0;
    Vector2 const s2 = p3 - p2;

    Vector2 const u = p0 - p2;

    float const ip = 1.f / (-s2.x * s1.y + s1.x * s2.y);

    float const s = (-s1.y * u.x + s1.x * u.y) * ip;
    float const t = ( s2.x * u.y - s2.y * u.x) * ip;

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        if (i) *i = p0 + (s1 * t);
        return true;
    }

    return false;
}
*/
func LineIntersection2(e1 Edge, e2 Edge) Vector {

    fmt.Printf("No no no! e1: %v, e2: %v\n", e1.Dir, e2.Dir)

    s1 := e1.Dir
    s2 := e2.Dir

    u  := Sub(e1.Pos, e2.Pos)

    ip := 1.0 / (-s2.X * s1.Y + s1.X * s2.Y)

    s  := (-s1.Y * u.X + s1.X * u.Y) * ip
    t  := ( s2.X * u.Y - s2.Y * u.X) * ip

    if s >= 0 && s <= 1 && t >= 0 && t <= 1 {
        return Add(e1.Pos, Mult(s1, t))
    }

    fmt.Println(s, t)

    fmt.Println("Lines do (still) NOT intersect!")
    return Vector{}
}

func LineIntersection3(e1 Edge, e2 Edge) Vector {

    path := geo.NewPath()
    path.Push(geo.NewPoint(float64(e1.Pos.X), float64(e1.Pos.Y)))
    path.Push(geo.NewPoint(float64(Add(e1.Pos, e1.Dir).X), float64(Add(e1.Pos, e1.Dir).Y)))

    line := geo.NewLine(geo.NewPoint(float64(e2.Pos.X), float64(e2.Pos.Y)), geo.NewPoint(float64(Add(e2.Pos, e2.Dir).X), float64(Add(e2.Pos, e2.Dir).Y)))

    // intersects does a simpler check for yes/no
    if path.Intersects(line) {
        fmt.Println("this shit intersects!")
        // intersection will return the actual points and places on intersection
        points, segments := path.Intersection(line)

        for i, _ := range points {
            fmt.Printf("Intersection %d at %v with path segment %d\n", i, points[i], segments[i][0])
            return Vector{float32(points[i][0]), float32(points[i][1]),0}
        }
    }

    fmt.Println("What the shit man????????????")
    return Vector{}
}

/*
func (l *Line) Interpolate(percent float64) *Point {
    return &Point{
        l.a[0] + percent*(l.b[0]-l.a[0]),
        l.a[1] + percent*(l.b[1]-l.a[1]),
    }
}
func LineIntersection4(e1 Edge, e2 Edge) Vector {
    den := (line.b[1]-line.a[1])*(l.b[0]-l.a[0]) - (line.b[0]-line.a[0])*(l.b[1]-l.a[1])
    U1 := (line.b[0]-line.a[0])*(l.a[1]-line.a[1]) - (line.b[1]-line.a[1])*(l.a[0]-line.a[0])
    U2 := (l.b[0]   -l.a[0])   *(l.a[1]-line.a[1]) - (l.b[1]   -l.a[1])   *(l.a[0]-line.a[0])

    if den == 0 {
        // collinear, all bets are off
        if U1 == 0 && U2 == 0 {
            return InfinityPoint
        }

        return nil
    }

    if U1/den < 0 || U1/den > 1 || U2/den < 0 || U2/den > 1 {
        return nil
    }

    return l.Interpolate(U1 / den)
}
*/
func LineIntersection4(e1 Edge, e2 Edge) (bool, Vector) {

    e1E := Add(e1.Pos, e1.Dir)
    e2E := Add(e2.Pos, e2.Dir)

    det := e1.Dir.Y * e2.Dir.X - e1.Dir.X * e2.Dir.Y
    s1 := (e1E.X-e1.Pos.X) * (e2.Pos.Y-e1.Pos.Y) - (e1E.Y-e1.Pos.Y) * (e2.Pos.X-e1.Pos.X)
    s2 := (e2E.X-e2.Pos.X) * (e2.Pos.Y-e1.Pos.Y) - (e2E.Y-e2.Pos.Y) * (e2.Pos.X-e1.Pos.X)

    if det <= 0.0000001 && det >= -0.0000001 {
        fmt.Println("It says - All bets are off. What are we doing now?")
        // collinear, all bets are off
        if s1 == 0 && s2 == 0 {
            fmt.Println(".... So maybe infinity? Hmpf")
            // Some kind of infinity stuff-Point??
            // I could just say, they intersect in like 1000-distance.
            return false, Vector{}
        }

        return false, Vector{}
    }

    if s1/det < 0 || s1/det > 1 || s2/det < 0 || s2/det > 1 {
        fmt.Println("I think, there is an intersection. But not within the edges given...")
        return false, Vector{}
    }

    return true, Add(e2.Pos, Mult(e2.Dir, s1/det))
}



