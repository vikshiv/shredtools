"""Vendored from mumemto.viz_mums — polygon geometry only (no matplotlib)."""


def points_to_poly(points):
    starts, ends = tuple(zip(*points))
    points = starts + ends[::-1]
    return points


def get_mum_polygons(mums, centering, color="#00A2FF", inv_color="red"):
    polygons = []
    colors = []
    for l, starts, strands in mums:
        inverted = not strands[0]
        points = []
        for idx, (x, strand) in enumerate(zip(starts, strands)):
            if x == -1:
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                continue
            points.append(((centering[idx] + x, idx), (centering[idx] + x + l, idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) >= 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors
