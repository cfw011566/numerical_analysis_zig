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
    std.mem.copyForwards(f64, Q_list, Q_input);
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
