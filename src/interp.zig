/// Interpolation
const std = @import("std");

/// Algorithm 3.1 Neville's Iterated Interpolation
/// x :
/// x_list :
/// Q_list :
pub fn neville(x: f64, x_list: []f64, Q_input: []f64, debug: bool) !f64 {
    const n = x_list.len;
    const allocator = std.heap.page_allocator;
    var Q_list = try allocator.alloc(f64, n);
    defer allocator.free(Q_list);
    @memcpy(Q_list, Q_input);
    for (1..n) |i| {
        if (debug) {
            std.debug.print("Q_list {any}\n", .{Q_list});
        }
        for (1..n - i + 1) |j| {
            const m = n - j;
            if (debug) {
                std.debug.print("m = {} n = {} i = {} j = {}\n", .{ m, n, i, j });
            }
            const x_j = x_list[m];
            const x_j_1 = x_list[m - i];
            const Q_i = Q_list[m];
            const Q_i_1 = Q_list[m - 1];
            const Q = ((x - x_j_1) * Q_i - (x - x_j) * Q_i_1) / (x_j - x_j_1);
            if (debug) {
                std.debug.print("Q = {d:.8}\n", .{Q});
            }
            Q_list[m] = Q;
        }
    }
    return Q_list[n - 1];
}

/// Algorithm 3.2 Newton's Divided-Difference Formula
/// x_list : x0, x1, ... xn
/// F_list : f(x0), f(x1), ... f(xn)
pub fn divided(x_list: []f64, F_input: []f64, F_output: []f64, debug: bool) !void {
    const n = x_list.len;
    @memcpy(F_output, F_input);
    for (1..n) |i| {
        if (debug) {
            std.debug.print("F {any}\n", .{F_output});
        }
        for (1..n - i + 1) |j| {
            const m = n - j;
            if (debug) {
                std.debug.print("m = {} n = {} i = {} j = {}\n", .{ m, n, i, j });
            }
            const x_j = x_list[m];
            const x_j_1 = x_list[m - i];
            const F_i = F_output[m];
            const F_i_1 = F_output[m - 1];
            const F = (F_i - F_i_1) / (x_j - x_j_1);
            if (debug) {
                std.debug.print("F = {d:.8}\n", .{F});
            }
            F_output[m] = F;
        }
    }
    if (debug) {
        std.debug.print("F_output {any}\n", .{F_output});
    }
}

pub fn divided_interp(x: f64, x_list: []f64, F_input: []f64, debug: bool) !f64 {
    const n = x_list.len;
    const allocator = std.heap.page_allocator;
    const F_output = try allocator.alloc(f64, n);
    defer allocator.free(F_output);
    try divided(x_list, F_input, F_output, debug);
    if (debug) {
        std.debug.print("F_output {any}\n", .{F_output});
    }
    var v: f64 = 0;
    for (0..n) |i| {
        var factor: f64 = 1.0;
        for (0..i) |j| {
            factor *= (x - x_list[j]);
            if (debug) {
                std.debug.print("i = {} j = {} x = {} factor = {d:.8}\n", .{ i, j, x_list[j], factor });
            }
        }
        v += factor * F_output[i];
    }
    return v;
}
