#include <iostream>
#include <math.h>
#include <random>
#include <thread>
#include <iomanip>

// OpenGL helpers (this should be put into one thing!)
#include "VertexBuffer/VertexBuffer.hpp"
#include "IndexBuffer/IndexBuffer.hpp"
#include "VertexArray/VertexArray.hpp"
#include "Shader/Shader.hpp"
#include "VertexBufferLayout/VertexBufferLayout.hpp"
#include "Rectangle/Rectangle.hpp"
#include "Circle/Circle.hpp"

// project-specific includes go here
#include "Graph/Graph.hpp"
#include "Utilities/Utilities.hpp"

SDL_Window* window;
SDL_GLContext gl_context;

const bool renderGraphics { true };

bool is_running;
bool paused { true };
bool step { false };

const float width { 768.0f };
const float height { 768.0f };

double DT { 0.025 };

// project-specific settings
unsigned int res { 12 };

std::random_device rd;

glm::fvec4 colorNode(float s)
{
    if (s > 0.0) { return glm::fvec4(0.0f, 1.0f, 1.0f, pow(abs(s),0.5)); }
    else { return glm::fvec4(1.0f, 0.0f, 0.0f, pow(abs(s),0.5)); }
}

void set_sdl_gl_attributes()
{
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
}

void evolveTime(SDL_Event& event)
{
    if (event.type == SDL_KEYDOWN)
    {
        if (event.key.keysym.scancode == SDL_SCANCODE_SPACE)
        {
            paused = paused ? false : true;
        }
        else if (event.key.keysym.scancode == SDL_SCANCODE_RIGHT)
        {
            step = true;
        }
        // else if (event.key.keysym.scancode == SDL_SCANCODE_LEFT)
        // {
        //     GLOBAL_TIME < DT ? GLOBAL_TIME = 0.0 : GLOBAL_TIME -= DT; // prevent negative time
        // }
    }
}

void printKey()
{
    std::cout << '\n' << "###################################" << '\n';
    std::cout << "KEY" << '\n';
    std::cout << "Spacebar" << '\t' << "Play/pause" << '\n';
    std::cout << "Right arrow" << '\t' << "Time forward" << '\n';
    // std::cout << "Left arrow" << '\t' << "Time reverse" << '\n';
    // std::cout << "Enter/Return" << '\t' << "Reset particles" << '\n';
    std::cout << "###################################" << '\n' << '\n';
}

int main()
{
    if (renderGraphics)
    {
        if(SDL_Init(SDL_INIT_EVERYTHING)==0)
        {
            std::cout<<"SDL2 initialized successfully."<<std::endl;
            set_sdl_gl_attributes();
            
            window = SDL_CreateWindow("Physarum Polycephalum", 0.0f, 0.0f, static_cast<int>(width), static_cast<int>(height), SDL_WINDOW_OPENGL);
            gl_context = SDL_GL_CreateContext(window);
            SDL_GL_SetSwapInterval(1);

            if(glewInit() == GLEW_OK)
            {
                std::cout << "GLEW initialization successful" << std::endl;
            }
            else
            {
                std::cout << "GLEW initialization failed" << std::endl;
                return -1;
            }

            std::cout << "Thread count: " << std::thread::hardware_concurrency() << '\n';

            GLCall(glEnable(GL_BLEND));
            GLCall(glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

            printKey();

            // buffer stuff and initialization happens here
            // make graph odd x odd so that there is always a central node
            Graph graph(rd(), width, height, 2 * res + 1);

            std::vector<Circle> circs;
            circs.reserve(graph.nodeCount());
            
            // create circles
            const std::vector<Node>& nodes { graph.nodes() };
            const Eigen::VectorXd& s { graph.getS() };
            for (unsigned int i = 0; i < nodes.size(); ++i)
            {
                // only create circle for sources/sinks
                if (s(i) != 0.0)
                {
                    const Node& node { nodes[i] };
                    float fs { static_cast<float>(s(i)) };
                    circs.emplace_back(Circle(node.pos, 0.2f * std::sqrt(width * height / static_cast<float>(graph.nodeCount())), colorNode(fs)));
                }
            }

            // prepare circles for graphics
            Circles circles(circs);

            VertexArray circ_VAO, line_VAO;

            // constructor automatically binds buffer
            VertexBuffer circ_VBO(circles.m_vertices.data(), static_cast<unsigned int>(circles.m_vertices.size() * sizeof(float)), GL_STATIC_DRAW);

            VertexBufferLayout circ_layout;
            for (int i = 0; i < 5; ++i)
            {
                circ_layout.push<float>(Circle::layout_descriptor[i]);
            }
            circ_VAO.addBuffer(circ_VBO, circ_layout);

            // constructor automatically binds buffer
            IndexBuffer circ_IBO(circles.m_indices.data(), static_cast<unsigned int>(circles.m_indices.size()));
            
            // const Eigen::SparseMatrix<double>& Laplacian { graph.getL() };
            const std::vector<Edge>& edges { graph.edges() };
            const Eigen::VectorXd& D { graph.getD() };
            const Eigen::VectorXd& Q { graph.getQ() };

            // TODO:
            // 1. add `line` class to OpenGL engine
            // 2. switch lines to quads that can change thickness (based off of D or Q)
            std::vector<float> lines;
            std::size_t edgeCount { graph.edgeCount() };
            lines.reserve(edgeCount);
            for (unsigned int i = 0; i < edgeCount; ++i)
            {
                const Edge& edge = edges[i];

                // node i
                lines.emplace_back(nodes[edge.i].pos[0]);
                lines.emplace_back(nodes[edge.i].pos[1]);

                lines.emplace_back(1.0f);
                lines.emplace_back(1.0f);
                lines.emplace_back(1.0f);
                lines.emplace_back(abs(edge.D));
                // OK to use `edge.D` here since it's synced with Qvec upon initialization

                // node j
                lines.emplace_back(nodes[edge.j].pos[0]);
                lines.emplace_back(nodes[edge.j].pos[1]);

                lines.emplace_back(1.0f);
                lines.emplace_back(1.0f);
                lines.emplace_back(1.0f);
                lines.emplace_back(abs(edge.D));
            }

            VertexBuffer line_VBO(lines.data(), static_cast<unsigned int>(lines.size() * sizeof(float)), GL_DYNAMIC_DRAW);

            VertexBufferLayout line_layout;
            line_layout.push<float>(2); // x,y
            line_layout.push<float>(4); // r,g,b,a
            line_VAO.addBuffer(line_VBO, line_layout);

            // set up MPV matrix
            glm::mat4 proj { glm::ortho(0.0f, width, 0.0f, height, -1.0f, 1.0f) };
            glm::mat4 view { glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0)) };
            glm::mat4 model { glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)) };
            glm::mat4 MVP { proj * view * model };

            // circ_shader stuff happens here
            Shader circ_shader("/Users/max/OpenGL_Framework/res/shaders", "circle_");
            circ_shader.bind();
            circ_shader.setUniformMatrix4fv("u_MVP", MVP);

            Shader line_shader("res/shaders", "colored_line_");
            line_shader.bind();
            line_shader.setUniformMatrix4fv("u_MVP", MVP);

            circ_VBO.unbind();
            circ_VAO.unbind();
            circ_IBO.unbind();
            circ_shader.unbind();

            line_VBO.unbind();
            line_VAO.unbind();
            line_shader.unbind();

            Renderer renderer;

            // Main loop
            SDL_Event event;
            is_running = true;

            while(is_running)
            {
                while(SDL_PollEvent(&event))
                {
                    if( event.type == SDL_QUIT || (event.type == SDL_KEYDOWN && event.key.keysym.scancode == SDL_SCANCODE_ESCAPE) )
                    {
                        is_running = false;
                    }

                    evolveTime(event);
                }

                // evolve time if unpaused
                if (not paused)
                {
                    graph.evolveGraph(DT);
                    // func wrapper for "write fitness to file"
                    // std::cout << "##############################################" << '\n';
                    // std::cout << "DE :" << '\t' << std::fixed << std::setprecision(9) << graph.dissipation() << '\n';
                    // std::cout << "TE :" << '\t' << graph.efficiency(D) << '\n';
                    if (graph.conductanceConverged() && graph.fitConverged())
                    {
                        // continue;
                        std::cout << "Conductance and fitness converged!" << '\n';
                        // paused = true;
                        // const std::string path { "/Users/max/TKN_Physarum/" };
                        // Utilities::exportCSV(path + "data", graph.sampleHSpec(1000, 1e-1));
                        // graph.printSpec(graph.sampleHSpec(50, 1e-2));
                        is_running = false;
                    }
                }
                if (step)
                {
                    // std::mt19937 rng(rd());
                    // std::cout << graph.probeHessian(rng, 1e-1) << '\n';
                    graph.evolveGraph(DT);
                    step = false;
                }

                // updating and rendering stuff happens here
                
                // update buffers (circle position, color, alpha)
                // * node values shouldn't change in current implementation *
                // circles.updateColors(values);
                // circ_VBO.updateBuffer(circles.m_vertices.data());

                // update line buffer
                for (unsigned int i = 0; i < edgeCount; ++i)
                {
                    // comment the two corresponding lines to have flow show up as that color
                    // red
                    // lines[12 * i + 2      ] = 1.0f - static_cast<float>(10 * abs(Q(i)));
                    // lines[12 * (i + 1) - 4] = 1.0f - static_cast<float>(10 * abs(Q(i)));
                    
                    // green
                    // lines[12 * i + 3      ] = 1.0f - static_cast<float>(10 * abs(Q(i)));
                    // lines[12 * (i + 1) - 3] = 1.0f - static_cast<float>(10 * abs(Q(i)));

                    // blue
                    lines[12 * i + 4      ] = 1.0f - static_cast<float>(10 * abs(Q(i)));
                    lines[12 * (i + 1) - 2] = 1.0f - static_cast<float>(10 * abs(Q(i)));

                    // alpha based on D
                    lines[12 * i + 5      ] = static_cast<float>(5 * pow(abs(D(i)), 0.5));
                    lines[12 * (i + 1) - 1] = static_cast<float>(5 * pow(abs(D(i)), 0.5));
                }
                line_VBO.updateBuffer(lines.data());

                renderer.clear();

                // lines first so that they "end" at the circle edge
                renderer.drawLines(line_VAO, line_VBO, line_shader);
                renderer.drawCircles(circ_VAO, circ_IBO, circ_shader);

                SDL_GL_SwapWindow(window);
            }

            SDL_Quit();
        }
        else
        {
            std::cout<<"SDL2 initialization failed."<<std::endl;
            return -1;
        }

        return 0;
    }
    else
    {
        for (int k = 0; k < 10; ++k)
        {
            // parallelize this!
            std::cout << "iteration " << k << '\n';
            Graph graph(rd(), width, height, 2 * res + 1);
            while (true)
            {
                graph.evolveGraph(DT);
                if (graph.conductanceConverged() && graph.fitConverged())
                {
                    std::cout << "Conductance and fitness converged!" << '\n';
                    const std::string path { "/Users/max/TKN_Physarum/kdata/" };
                    Utilities::exportCSV(path + std::to_string(k), graph.sampleHSpec(1, 1e-2));
                    break;
                }
            }
    }
    }
}
