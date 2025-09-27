#include "pr.hpp"
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

inline float asymmetric_gaussian(float x, float mean, float sigma, float a) {
    float k = (x - mean) / (sigma + a * (x - mean));
    return expf(-k * k);
}
inline float cmf_x_CIE_2006_10deg(float x) {
    float a = 0.42f * asymmetric_gaussian(x, 445.5793762207031f, 31.10230255126953f, 0.06199074536561966f);
    float b = 1.16f * asymmetric_gaussian(x, 594.541748046875f, 48.5798225402832f, -0.04986543208360672f);
    return a + b - a * b * 42.51573944091797f;
}
inline float cmf_y_CIE_2006_10deg(float x) {
    return 1.0f * asymmetric_gaussian(x, 556.8382568359375f, 66.54190826416016f, -0.026490915566682816f);
}
inline float cmf_z_CIE_2006_10deg(float x) {
    return 2.146832f * asymmetric_gaussian(x, 445.9251708984375f, 30.91781997680664f, 0.08379142731428146f);
}

// importance sampling
inline float logistic_pdf(float x, float s)
{
    float k = expf(-fabsf(x) / s);
    return s * k / ((1.0 + k) * (1.0 + k));
}
inline float logistic_cdf(float x, float s)
{
    return 1.0f / (1.0f + expf(-x / s));
}
inline float inverse_logistic_cdf(float u, float s)
{
    if (0.99999994f < u) { u = 0.99999994f; }
    if (u < 1.175494351e-38f) { u = 1.175494351e-38f; }
    return -s * logf(1.0f / u - 1.0f);
}
inline float trimmed_logistic_pdf(float x, float s, float a, float b)
{
    return logistic_pdf(x, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
}

inline float cmf_y_pdf_CIE_2006_10deg(float x) {
    float sx = x - 554.270751953125f;
    float s = 26.879621505737305f;
    float a = -164.270751953125;
    float b = 275.729248046875;
    return logistic_pdf(sx, 26.879621505737305f) / (logistic_cdf(b, s) - logistic_cdf(a, s));
}

inline float cmf_y_sample_CIE_2006_10deg(float u) {
    float s = 26.879621505737305f;
    float a = -164.270751953125;
    float b = 275.729248046875;
    float Pa = logistic_cdf(a, s);
    float Pb = logistic_cdf(b, s);
    return inverse_logistic_cdf(Pa + (Pb - Pa) * u, s) + 554.270751953125f;
}

int main() {
    using namespace pr;

    std::vector<float> wavelength;
    std::vector<float> Xs;
    std::vector<float> Ys;
    std::vector<float> Zs;

    std::ifstream csv("../lin2012xyz10e_1_7sf.csv");
    std::string line;
    while (std::getline(csv, line)) {
        std::stringstream ss(line);

        float wave, x, y, z;
        sscanf(line.c_str(), "%f,%f,%f,%f", &wave, &x, &y, &z);
        wavelength.push_back(wave);
        Xs.push_back(x);
        Ys.push_back(y);
        Zs.push_back(z);
    }

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 5, 0, 10 };
    camera.lookat = { 5, 0, 0 };

    double e = GetElapsedTime();

    bool showOriginalData = false;
    bool showX = true;
    bool showY = true;
    bool showZ = true;

    bool showPDF = false;
    bool showSampledHistogram = false;

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        ClearBackground(0.1f, 0.1f, 0.1f, 1);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XY, 1.0f, 20, { 128, 128, 128 });
        DrawXYZAxis(1.0f);

        for (int nm = 0; nm <= 800; nm += 100)
        {
            float x = nm / 100.0f;
            DrawText({ x, 0.0f, 0.0f }, std::to_string(nm) + " nm", 12);
        }

        // Plot
        if (showOriginalData)
        {
            for (int i = 0; i < wavelength.size(); i++)
            {
                float nm = wavelength[i];
                float x = nm / 100.0f;
                DrawPoint({ x, Xs[i], 0.0f }, { 128, 0, 0 }, 4);
                DrawPoint({ x, Ys[i], 0.0f }, { 0, 128, 0 }, 4);
                DrawPoint({ x, Zs[i], 0.0f }, { 0, 0, 128 }, 4);
            }
        }

        // Fitted
        int N = 1000;

        if (showX)
        {
            PrimBegin(PrimitiveMode::LineStrip, 1);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_x_CIE_2006_10deg(nm), 0 }, { 255, 0, 0 });
            }
            PrimEnd();
        }

        if (showY)
        {
            PrimBegin(PrimitiveMode::LineStrip, 1);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_y_CIE_2006_10deg(nm), 0 }, { 0, 255, 0 });
            }
            PrimEnd();
        }

        if (showZ)
        {
            PrimBegin(PrimitiveMode::LineStrip, 1);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_z_CIE_2006_10deg(nm), 0 }, { 0, 0, 255 });
            }
            PrimEnd();
        }

        // PDF
        static float pdf_view_scale = 0.2f;

        if (showPDF)
        {
            PrimBegin(PrimitiveMode::LineStrip, 1);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_y_pdf_CIE_2006_10deg(nm) * pdf_view_scale, 0 }, { 255, 255, 255 });
            }
            PrimEnd();
        }

        // sampling
        static std::vector<int> buckets(900);
        static int maxCount = 0;
        if (showSampledHistogram)
        {
            PCG rng;

            if (maxCount == 0)
            {
                for (int i = 0; i < 100000000; i++)
                {
                    float sampled = cmf_y_sample_CIE_2006_10deg(rng.uniformf());
                    int index = (int)roundf(sampled);
                    int cur = ++buckets[index];
                    maxCount = std::max(maxCount, cur);
                }
            }

            PrimBegin(PrimitiveMode::LineStrip, 2);
            for (int nm = 0; nm < buckets.size(); nm++)
            {
                float x = nm / 100.0f;
                float y = (float)buckets[nm] / maxCount;
                PrimVertex({ x, y, 0 }, { 255, 255, 0 });
            }
            PrimEnd();
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        ImGui::Checkbox("show original data", &showOriginalData);

        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(255, 0, 0, 255));
        ImGui::Checkbox("showX", &showX);
        ImGui::PopStyleColor();

        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(0, 255, 0, 255));
        ImGui::Checkbox("showY", &showY);
        ImGui::PopStyleColor();

        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(0, 0, 255, 255));
        ImGui::Checkbox("showZ", &showZ);
        ImGui::PopStyleColor();

        ImGui::Checkbox("show PDF(y)", &showPDF);
        if (showPDF)
        {
            ImGui::SliderFloat("pdf view scale", &pdf_view_scale, 0, 1);
        }
        ImGui::Checkbox("show Sampled Histogram(y)", &showSampledHistogram);

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
