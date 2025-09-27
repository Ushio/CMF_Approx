#include "pr.hpp"
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

namespace CIE_1931_2deg
{
    inline float asymmetric_gaussian(float x, float mean, float sigma, float a) {
        float k = (x - mean) / (sigma + a * (x - mean));
        return expf(-k * k);
    }
    inline float cmf_x(float x) {
        float a = 0.42f * asymmetric_gaussian(x, 447.95562744140625f, 38.40695571899414f, 0.40412160754203796f);
        float b = 1.64f * asymmetric_gaussian(x, 588.0320434570312f, 55.6748046875f, -0.1073513776063919f);
        return a + b - a * b * 8.491180419921875f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.5616455078125f, 59.5950927734375f, 0.056370146572589874f);
    }
    inline float cmf_z(float x) {
        return 1.7829682f * asymmetric_gaussian(x, 447.90704345703125f, 32.452816009521484f, 0.12668778002262115f);
    }

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

    inline float cmf_y_pdf(float x) {
        float sx = x - 559.8692016601562f;
        float s = 23.981721878051758f;
        float a = -169.86920166015625;
        float b = 270.13079833984375;
        return logistic_pdf(sx, 23.981721878051758f) / (logistic_cdf(b, s) - logistic_cdf(a, s));
    }

    inline float cmf_y_sample(float u) {
        float s = 23.981721878051758f;
        float a = -169.86920166015625;
        float b = 270.13079833984375;
        float Pa = logistic_cdf(a, s);
        float Pb = logistic_cdf(b, s);
        return inverse_logistic_cdf(Pa + (Pb - Pa) * u, s) + 559.8692016601562f;
    }
}

namespace CIE_2015_10deg
{
    inline float asymmetric_gaussian(float x, float mean, float sigma, float a) {
        float k = (x - mean) / (sigma + a * (x - mean));
        return expf(-k * k);
    }
    inline float cmf_x(float x) {
        float a = 0.42f * asymmetric_gaussian(x, 445.5875549316406f, 31.161394119262695f, 0.065997414290905f);
        float b = 1.16f * asymmetric_gaussian(x, 594.5599975585938f, 48.59465408325195f, -0.04558335989713669f);
        return a + b - a * b * 42.55435562133789f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.8383178710938f, 66.54190826416016f, -0.026492968201637268f);
    }
    inline float cmf_z(float x) {
        return 2.146832f * asymmetric_gaussian(x, 445.9251708984375f, 30.91781997680664f, 0.08379140496253967f);
    }

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

    inline float cmf_y_pdf(float x) {
        float sx = x - 554.270751953125f;
        float s = 26.879621505737305f;
        float a = -164.270751953125;
        float b = 275.729248046875;
        return logistic_pdf(sx, 26.879621505737305f) / (logistic_cdf(b, s) - logistic_cdf(a, s));
    }

    inline float cmf_y_sample(float u) {
        float s = 26.879621505737305f;
        float a = -164.270751953125;
        float b = 275.729248046875;
        float Pa = logistic_cdf(a, s);
        float Pb = logistic_cdf(b, s);
        return inverse_logistic_cdf(Pa + (Pb - Pa) * u, s) + 554.270751953125f;
    }
}

int main() {
    using namespace pr;

    // 1931
    //using namespace CIE_1931_2deg;
    //const char* dataName = "../cie_2_1931.csv";

    // 2015
    using namespace CIE_2015_10deg;
    const char* dataName = "../CIE 2015 10 Degree Standard Observer.csv";


    // loading
    std::vector<float> wavelength;
    std::vector<float> Xs;
    std::vector<float> Ys;
    std::vector<float> Zs;

    std::ifstream csv(dataName);
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
    camera.origin = { 5, 1, 4 };
    camera.lookat = { 5, 1, 0 };

    double e = GetElapsedTime();

    bool showOriginalData = true;
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
            PrimBegin(PrimitiveMode::LineStrip, 3);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_x(nm), 0 }, { 255, 0, 0 });
            }
            PrimEnd();
        }

        if (showY)
        {
            PrimBegin(PrimitiveMode::LineStrip, 3);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_y(nm), 0 }, { 0, 255, 0 });
            }
            PrimEnd();
        }

        if (showZ)
        {
            PrimBegin(PrimitiveMode::LineStrip, 3);
            for (int i = 0; i < N; i++)
            {
                float nm = glm::mix(390.0f, 830.0f, (float)i / N);
                float x = nm / 100.0f;
                PrimVertex({ x, cmf_z(nm), 0 }, { 0, 0, 255 });
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
                PrimVertex({ x, cmf_y_pdf(nm) * pdf_view_scale, 0 }, { 255, 255, 255 });
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
                    float sampled = cmf_y_sample(rng.uniformf());
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

        ImGui::SetNextWindowSize({ 400, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());


        ImGui::Text("%s", dataName);
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
