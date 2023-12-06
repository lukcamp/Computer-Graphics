#pragma once

#include <Eigen/Core>

class VertexAttributes
{
public:
    VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1)
    {
        position << x, y, z, w;
        normal << 0, 0, 0, 0;
        diffuse_color << 0, 0, 0, 1; // RGBA color
        specular_color << 0, 0, 0, 1; // RGBA color
        specular_exponent = 0;
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes &a,
        const VertexAttributes &b,
        const VertexAttributes &c,
        const double alpha,
        const double beta,
        const double gamma)
    {
        VertexAttributes r;
        r.position = alpha * a.position + beta * b.position + gamma * c.position;
        return r;
    }

    Eigen::Vector4d position;
    Eigen::Vector4d normal;
    Eigen::Vector3d color;
    Eigen::Vector4d diffuse_color;
    Eigen::Vector4d specular_color;
    double specular_exponent;
};

class FragmentAttributes
{
public:
    FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1, double d = 0)
    {
        color << r, g, b, a;
        depth = d;
    }

    Eigen::Vector4d color;
    double depth;
};

class FrameBufferAttributes
{
public:
    FrameBufferAttributes(double r = 0, double g = 0, double b = 0, double a = 255, double d = 0)
    {
        color << r, g, b, a;
        depth = d;
    }

    Eigen::Matrix<double, 4, 1> color;
    double depth;
};

class UniformAttributes
{
public:
    Eigen::Matrix4d view_matrix;
    Eigen::Matrix4d projection_matrix;
};
