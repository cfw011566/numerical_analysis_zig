/// Root finding
const std = @import("std");

const RootError = error{ WrongEndPoints, OverIteration, NotInSection, NoRoot };
const Function = *const fn (x: f64) f64;
const Function2 = *const fn (x: std.math.Complex(f64)) std.math.Complex(f64);

/// Algorithm 2.1 Bisection
/// a0, b0 : end points
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn bisection(f: Function, a0: f64, b0: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    if (a0 >= b0) {
        return error.WrongEndPoints;
    }
    var a = a0;
    var b = b0;
    var fa = f(a);
    var fb = f(b);
    if (std.math.sign(fa) * std.math.sign(fb) > 0) {
        return error.NotInSection;
    }

    var n: usize = 1;
    while (n <= N0) : (n += 1) {
        const c0 = (b - a) / 2;
        // p.52 of burden
        // const c = (a + b) / 2;
        const c = a + c0;
        const relative_error = c0 / c;
        const fc = f(c);
        if (debug) {
            std.debug.print("{d:>3} {d:10.6} {d:10.6} {d:10.6} {d:10.6} {d:10.6}\n", .{ n, a, b, c, fc, relative_error });
        }
        if (fc == 0.0 or c0 < TOL) {
            return c;
        }
        if (std.math.sign(fc) * std.math.sign(fa) < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return error.OverIteration;
}

/// p0 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn fixed_point(f: Function, p0: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var p = p0;
    var n: usize = 0;
    var p1 = p;
    while (n <= N0) : (n += 1) {
        p1 = f(p);
        if (debug) {
            std.debug.print("{d:>3} {d:.10} {d:.10}\n", .{ n, p, p1 - p });
        }
        if (std.math.approxEqAbs(f64, p, p1, TOL)) {
            return p1;
        }
        p = p1;
    }
    return error.OverIteration;
}

/// p0 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn newton(f: Function, f1: Function, p0: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var p = p0;
    var p1 = p;
    var n: usize = 0;
    while (n <= N0) : (n += 1) {
        p1 = p - f(p) / f1(p);
        if (debug) {
            std.debug.print("{d:>3} {d:.10} {d:.10}\n", .{ n, p, p1 - p });
        }
        if (std.math.approxEqAbs(f64, p, p1, TOL)) {
            return p1;
        }
        p = p1;
    }
    return error.OverIteration;
}

/// p0, p1 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn secant(f: Function, p0: f64, p1: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var a = p0;
    var b = p1;
    var fa = f(a);
    var fb = f(b);
    var n: usize = 2;
    while (n <= N0) : (n += 1) {
        const c = b - fb * (b - a) / (fb - fa);
        if (debug) {
            std.debug.print("{d:>3} {d:.10} {d:.10}\n", .{ n, c, c - b });
        }
        if (std.math.approxEqAbs(f64, c, b, TOL)) {
            return c;
        }
        a = b;
        fa = fb;
        b = c;
        fb = f(c);
    }
    return error.OverIteration;
}

/// p0, p1 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn false_position(f: Function, p0: f64, p1: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var a = p0;
    var b = p1;
    var fa = f(a);
    var fb = f(b);
    var i: usize = 2;
    while (i <= N0) : (i += 1) {
        const c = b - fb * (b - a) / (fb - fa);
        if (debug) {
            std.debug.print("{d:>3} {d:.10} {d:.10}\n", .{ i, c, c - b });
        }
        if (std.math.approxEqAbs(f64, c, b, TOL)) {
            return c;
        }
        const fc = f(c);
        if (std.math.sign(fc) * std.math.sign(fb) < 0) {
            a = b;
            fa = fb;
        }
        b = c;
        fb = fc;
    }
    return error.OverIteration;
}

/// p0 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn steffensen(g: Function, p: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var p0 = p;
    var n: usize = 0;
    while (n <= N0) : (n += 1) {
        const p1 = g(p0);
        const p2 = g(p1);
        const pp = p0 - (p1 - p0) * (p1 - p0) / (p2 - p1 - p1 + p0);
        if (debug) {
            std.debug.print("{d:>3} {d:.10} {d:.10} {d:.10} {d:.10}\n", .{ n, p0, p1, p2, pp });
        }
        if (std.math.approxEqAbs(f64, pp, p0, TOL)) {
            return pp;
        }
        p0 = pp;
    }
    return error.OverIteration;
}

/// p0, p1, p2 : initial approximation
/// TOL : tolerance
/// N0 : maximum number of iterations
pub fn muller(f: Function2, p_0: f64, p_1: f64, p_2: f64, TOL: f64, N0: usize, debug: bool) !f64 {
    var p0 = std.math.Complex(f64).init(p_0, 0);
    var p1 = std.math.Complex(f64).init(p_1, 0);
    var p2 = std.math.Complex(f64).init(p_2, 0);
    var h1 = p1.sub(p0);
    var h2 = p2.sub(p1);
    var delta1 = f(p1).sub(f(p0)).div(h1);
    var delta2 = f(p2).sub(f(p1)).div(h2);
    var d = delta2.sub(delta1).div(h2.add(h1));
    var n: usize = 3;
    const four = std.math.Complex(f64).init(4.0, 0.0);
    const minus_two = std.math.Complex(f64).init(-2.0, 0.0);
    while (n <= N0) : (n += 1) {
        var b = delta2.add(h2.mul(d));
        const D_2 = b.mul(b).sub((f(p2).mul(four)).mul(d));
        const D = std.math.complex.sqrt(D_2);
        const b_minus_D = b.sub(D);
        const b_plus_D = b.add(D);
        const E = if (std.math.complex.abs(b_minus_D) < std.math.complex.abs(b_plus_D)) b.add(D) else b.sub(D);
        const h = minus_two.mul(f(p2)).div(E);
        const p = p2.add(h);
        if (debug) {
            std.debug.print("{d:>3} p = {} f(p) = {}\n", .{ n, p, f(p) });
        }
        if (std.math.complex.abs(h) < TOL) {
            if (p.im == 0.0) {
                return p.re;
            } else {
                return error.NoRoot;
            }
        }
        p0 = p1;
        p1 = p2;
        p2 = p;
        h1 = p1.sub(p0);
        h2 = p2.sub(p1);
        delta1 = f(p1).sub(f(p0)).div(h1);
        delta2 = f(p2).sub(f(p1)).div(h2);
        d = delta2.sub(delta1).div(h2.add(h1));
    }
    return error.OverIteration;
}
