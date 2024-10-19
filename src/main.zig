const std = @import("std");

const root = @import("root.zig");
const interp = @import("interp.zig");

pub fn main() !void {
    // try chapter2();
    try chapter3();
}

fn chapter2() !void {
    const r = try root.bisection(fn_ch_2_1_example_1, 1, 2, 1e-4, 20, false);
    std.debug.print("chapter 2.1 example 1 = {d:.6}\n", .{r});
    const r2 = try root.fixed_point(fn_ch_2_2, 1.5, 1e-8, 20, false);
    std.debug.print("chapter 2.2 = {d:.8}\n", .{r2});
    const p0: f64 = std.math.pi / 4.0;
    std.debug.print("chapter 2.3 example 1 (fixed point) = {d:.10}\n", .{try root.fixed_point(fn_ch_2_3_ex_1_fp, p0, 1e-4, 20, true)});
    std.debug.print("chapter 2.3 example 1 (newton) = {d:.10}\n", .{try root.newton(fn_ch_2_3_ex_1, fn_ch_2_3_ex_1_p, p0, 1e-8, 20, true)});
    std.debug.print("chapter 2.3 example 1 (secant) = {d:.10}\n", .{try root.secant(fn_ch_2_3_ex_1, 0.5, p0, 1e-8, 20, true)});
    std.debug.print("chapter 2.3 example 1 (false position) = {d:.10}\n", .{try root.false_position(fn_ch_2_3_ex_1, 0.5, p0, 1e-8, 20, true)});
    std.debug.print("chapter 2.4 example 1 (newton) = {d:.10}\n", .{try root.newton(fn_ch_2_4_ex_1, fn_ch_2_4_ex_1_p, 1.0, 1e-6, 20, true)});
    std.debug.print("chapter 2.5 illustration 1 (steffensen) = {d:.10}\n", .{try root.steffensen(fn_ch_2_5_il_1, 1.5, 1e-8, 20, true)});

    //    std.debug.print("chapter 2.6 illustration (muller) = {d:.10}\n", .{try root.muller(fn_ch_2_6_il, 0.5, -0.5, 0, 1e-8, 20, false)});
    std.debug.print("chapter 2.6 illustration (muller) = {d:.10}\n", .{try root.muller(fn_ch_2_6_il, 0.5, 1.0, 1.5, 1e-8, 20, false)});
    std.debug.print("chapter 2.6 illustration (muller) = {d:.10}\n", .{try root.muller(fn_ch_2_6_il, 1.5, 2.0, 2.5, 1e-8, 20, false)});
}

fn fn_ch_2_1_example_1(x: f64) f64 {
    return x * x * x + 4 * x * x - 10;
}

fn fn_ch_2_2(x: f64) f64 {
    return x - (x * x * x + 4 * x * x - 10) / (3 * x * x + 8 * x);
}

fn fn_ch_2_3_ex_1_fp(x: f64) f64 {
    return std.math.cos(x);
}

fn fn_ch_2_3_ex_1(x: f64) f64 {
    return std.math.cos(x) - x;
}

fn fn_ch_2_3_ex_1_p(x: f64) f64 {
    return -std.math.cos(x) - 1;
}

fn fn_ch_2_4_ex_1(x: f64) f64 {
    return std.math.exp(x) - x - 1;
}

fn fn_ch_2_4_ex_1_p(x: f64) f64 {
    return std.math.exp(x) - 1;
}

fn fn_ch_2_5_il_1(x: f64) f64 {
    return std.math.sqrt(10.0 / (x + 4.0));
}

fn fn_ch_2_6_il(x: std.math.Complex(f64)) std.math.Complex(f64) {
    var three = std.math.Complex(f64).init(3, 0);
    var one = std.math.Complex(f64).init(1, 0);
    var square = x;
    square = square.mul(square);
    var quad = square;
    quad = quad.mul(quad);
    var trip = square;
    trip = trip.mul(x);
    return one.add(x).add(square).sub(three.mul(trip)).add(quad);
}

fn chapter3() !void {
    try chapter3_1();
    try chapter3_2();
    try chapter3_exercise1();
    try chapter3_3_2_example1();
}

fn chapter3_1() !void {
    const x_input = [_]f64{ 1.0, 1.3, 1.6, 1.9, 2.2 };
    const Q_input = [_]f64{ 0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623 };
    const Q = try interp.neville(1.5, @constCast(x_input[0..]), @constCast(Q_input[0..]), false);
    std.debug.print("Q = {d:.7}\n", .{Q});
}

fn chapter3_2() !void {
    const x_input = [_]f64{ 2.0, 2.2, 2.3 };
    const Q_input = [_]f64{ 0.6931, 0.7885, 0.8329 };
    const Q = try interp.neville(2.1, @constCast(x_input[0..]), @constCast(Q_input[0..]), false);
    std.debug.print("Q = {d:.4}\n", .{Q});
}

fn chapter3_exercise1() !void {
    const x = 8.4;
    const x_input = [_]f64{ 8.1, 8.3, 8.6, 8.7 };
    const Q_input = [_]f64{ 16.94410, 17.56492, 18.50515, 18.82091 };
    const Q = try interp.neville(x, @constCast(x_input[0..]), @constCast(Q_input[0..]), false);
    std.debug.print("Q = {d:.7}\n", .{Q});
    const v = try interp.divided_interp(x, @constCast(x_input[0..]), @constCast(Q_input[0..]), false);
    std.debug.print("v = {d:.7}\n", .{v});
}

fn chapter3_3_2_example1() !void {
    const x_input = [_]f64{ 1.0, 1.3, 1.6, 1.9, 2.2 };
    const F_input = [_]f64{ 0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623 };
    var v = try interp.divided_interp(1.5, @constCast(x_input[0..]), @constCast(F_input[0..]), false);
    std.debug.print("f(1.5) = {d:.7}\n", .{v});
    v = try interp.divided_interp(1.1, @constCast(x_input[0..]), @constCast(F_input[0..]), false);
    std.debug.print("f(1.1) = {d:.7}\n", .{v});
    v = try interp.divided_interp(2.0, @constCast(x_input[0..]), @constCast(F_input[0..]), false);
    std.debug.print("f(2.0) = {d:.7}\n", .{v});
}
