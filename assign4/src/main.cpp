// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5;       //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    //TODO: setup uniform

    //TODO: setup camera, compute w, u, v

    Vector3d w = -camera_gaze.normalized();
    Vector3d u = camera_top.cross(w).normalized();
    Vector3d v = w.cross(u);

    //TODO: compute the camera transformation

    Matrix4d view_matrix;
    view_matrix <<
        u.x(), u.y(), u.z(), -u.dot(camera_position),
        v.x(), v.y(), v.z(), -v.dot(camera_position),
        w.x(), w.y(), w.z(), -w.dot(camera_position),
        0,     0,     0,     1;

    //TODO: setup projection matrix

    Matrix4d P;
    if (is_perspective)
    {
        double t = near_plane * tan(field_of_view / 2.0);
        double r = t * aspect_ratio;
        P <<
            near_plane / r, 0, 0, 0,
            0, near_plane / t, 0, 0,
            0, 0, -(far_plane + near_plane) / (far_plane - near_plane), -(2 * far_plane * near_plane) / (far_plane - near_plane),
            0, 0, -1, 0;
    }
    else
    {
        double t = 1.5;
        double l = -1.5;
        double r = 1.5;
        double b = -1.5;
        P <<
            2 / (r - l), 0, 0, -(r + l) / (r - l),
            0, 2 / (t - b), 0, -(t + b) / (t - b),
            0, 0, -2 / (far_plane - near_plane), -(far_plane + near_plane) / (far_plane - near_plane),
            0, 0, 0, 1;
    }

    uniform.view_matrix = view_matrix;
    uniform.projection_matrix = P;
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        VertexAttributes transformed_va;
        transformed_va.position = (uniform.projection_matrix * uniform.view_matrix * va.position.homogeneous()).hnormalized();
      // TODO: Pass any other necessary attributes to the fragment shader
      // For simple rendering, you might not need additional attributes.
        return transformed_va;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: build the vertex attributes from vertices and facets


    for (int i = 0; i < facets.rows(); ++i) {
        for (int j = 0; j < 3; ++j) { // Iterate over the three vertices of each triangle
            VertexAttributes va;
            int v_idx = facets(i, j);
            va.position = vertices.row(v_idx);
            // TODO: Pass other attributes like normals, colors, etc. if needed
            vertex_attributes.push_back(va);
        }
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

Matrix4d compute_rotation(const double alpha)
{
    //TODO: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    Matrix4d res = Matrix4d::Identity();

    double cos_alpha = cos(alpha);
    double sin_alpha = sin(alpha);

    res(0, 0) = cos_alpha;
    res(0, 2) = sin_alpha;
    res(2, 0) = -sin_alpha;
    res(2, 2) = cos_alpha;

    return res;
}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    Matrix4d trafo = compute_rotation(alpha);

    program.VertexShader = [&trafo](const VertexAttributes &va, const UniformAttributes &uniform) {

    //program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader

        VertexAttributes transformed_va;
        transformed_va.position = (trafo * va.position.homogeneous()).hnormalized();

        return transformed_va;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;

    //TODO: generate the vertex attributes for the edges and rasterize the lines
    //TODO: use the transformation matrix

    for (int i = 0; i < facets.rows(); ++i) {
            for (int j = 0; j < 3; ++j) { //iterate over the three vertices of each triangle
                VertexAttributes va;
                int v_idx = facets(i, j);
                va.position = vertices.row(v_idx);
                vertex_attributes.push_back(va);
            }
      }

    /*for (int i = 0; i < facets.rows(); ++i) {
        int v0_idx = facets(i, 0);
        int v1_idx = facets(i, 1);
        int v2_idx = facets(i, 2);

        VertexAttributes v0;
        v0.position = vertices.row(v0_idx);
        vertex_attributes.push_back(v0);

        VertexAttributes v1;
        v1.position = vertices.row(v1_idx);
        vertex_attributes.push_back(v1);

        VertexAttributes v2;
        v2.position = vertices.row(v2_idx);
        vertex_attributes.push_back(v2);
    }*/

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}

void get_shading_program(Program &program)
{
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: transform the position and the normal
        //TODO: compute the correct lighting
        VertexAttributes transformed_va;

        transformed_va.position = (uniform.projection_matrix * uniform.view_matrix * va.position.homogeneous()).hnormalized();
        transformed_va.normal = (uniform.view_matrix * va.normal).normalized();

        Vector3d light_direction = (uniform.view_matrix * (light_positions[0] - va.position.head<3>())).normalized();
        double diffuse_intensity = std::max(0.0, va.normal.dot(light_direction));
        Vector3d view_direction = (uniform.view_matrix * -va.position.head<3>()).normalized();
        Vector3d reflection_direction = (-light_direction) - 2.0 * (-light_direction).dot(va.normal) * va.normal;
        double specular_intensity = pow(std::max(0.0, reflection_direction.dot(view_direction)), va.specular_exponent);

        // Compute final vertex color based on lighting
        Vector3d diffuse_color = va.diffuse_color * diffuse_intensity * light_colors[0];
        Vector3d specular_color = va.specular_color * specular_intensity * light_colors[0];

        transformed_va.color = diffuse_color + specular_color;

        return transformed_va;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: create the correct fragment
        return FragmentAttributes(va.color[0], va.color[1], va.color[2]);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: implement the depth check
        if (fa.depth < previous.depth){
        // If the new fragment is closer, use its color
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], 1.0, fa.depth);
        }
        else return previous;

        //return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);
    Eigen::Matrix4d trafo = compute_rotation(alpha);

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: compute the normals


    for (int i = 0; i < facets.rows(); ++i) {
        const Eigen::Vector3d &v0 = vertices.row(facets(i, 0));
        const Eigen::Vector3d &v1 = vertices.row(facets(i, 1));
        const Eigen::Vector3d &v2 = vertices.row(facets(i, 2));
        Eigen::Vector3d face_normal = (v1 - v0).cross(v2 - v0).normalized();

        VertexAttributes va0;
        va0.position = v0;
        va0.normal = face_normal;
        vertex_attributes.push_back(va0);

        VertexAttributes va1;
        va1.position = v1;
        va1.normal = face_normal;
        vertex_attributes.push_back(va1);

        VertexAttributes va2;
        va2.position = v2;
        va2.normal = face_normal;
        vertex_attributes.push_back(va2);
    }

    //TODO: set material colors
    for (VertexAttributes &va : vertex_attributes) {
        va.diffuse_color = obj_diffuse_color;
        va.specular_color = obj_specular_color;
        va.specular_exponent = obj_specular_exponent;
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: compute the vertex normals as vertex normal average
    std::vector<Eigen::Vector3d> vertex_normals(vertices.rows(), Eigen::Vector3d::Zero());

    for (int i = 0; i < facets.rows(); ++i) {
        const Eigen::Vector3d &v0 = vertices.row(facets(i, 0));
        const Eigen::Vector3d &v1 = vertices.row(facets(i, 1));
        const Eigen::Vector3d &v2 = vertices.row(facets(i, 2));
        Eigen::Vector3d face_normal = (v1 - v0).cross(v2 - v0);

        vertex_normals[facets(i, 0)] += face_normal;
        vertex_normals[facets(i, 1)] += face_normal;
        vertex_normals[facets(i, 2)] += face_normal;
    }

    for (int i = 0; i < vertex_normals.size(); ++i) {
        vertex_normals[i].normalize();
    }

    for (int i = 0; i < vertices.rows(); ++i) {
        VertexAttributes va;
        va.position = vertices.row(i);
        va.normal = Eigen::Vector4d(vertex_normals[i][0], vertex_normals[i][1], vertex_normals[i][2], 0.0); // Convert to homogeneous coordinates
        va.diffuse_color = obj_diffuse_color;
        va.specular_color = obj_specular_color;
        va.specular_exponent = obj_specular_exponent;
        vertex_attributes.push_back(va);
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    //TODO: add the animation

    return 0;
}
