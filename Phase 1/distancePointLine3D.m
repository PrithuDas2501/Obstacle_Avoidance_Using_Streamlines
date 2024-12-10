function distance = distancePointLine3D(point, linePoint, lineDirection)
    % point: [x0, y0, z0]
    % linePoint: [x1, y1, z1] (point on the line)
    % lineDirection: [a, b, c] (direction vector of the line)

    % Calculate vector from line point to the given point
    vectorToLine = point - linePoint;

    % Calculate cross product
    crossProduct = cross(vectorToLine, lineDirection);

    % Calculate distance
    distance = norm(crossProduct) / norm(lineDirection);
end