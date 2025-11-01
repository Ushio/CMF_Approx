#include "pr.hpp"
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

namespace CIE_1931_2deg
{
    inline float asymmetric_gaussian(float x, float mean, float sigma, float a) {
        float denom = sigma + a * (x - mean);
        if (denom < 1.0e-15f) { denom = 1.0e-15f; }
        if (sigma * 2.0f < denom) { denom = sigma * 2.0f; }
        float k = (x - mean) / denom;
        return expf(-k * k);
    }
    inline float cmf_x(float x) {
        float a = 0.37f * asymmetric_gaussian(x, 445.8890380859375f, 32.71352767944336f, 0.2403123378753662f);
        float b = 1.113f * asymmetric_gaussian(x, 593.9199829101562f, 51.980140686035156f, -0.06552795320749283f);
        return (a + b - a * b * 21.016616821289062f) * 1.0057787153642084f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.5616455078125f, 59.5950927734375f, 0.056370146572589874f) * 1.0067831838201629f;
    }
    inline float cmf_z(float x) {
        return 1.7829682f * asymmetric_gaussian(x, 447.90643310546875f, 32.452659606933594f, 0.12635648250579834f) * 1.0171607842044876f;
    }

    inline float logistic_pdf(float x, float s)
    {
        float k = expf(-fabsf(x) / s);
        return k / ((1.0 + k) * (1.0 + k) * s);
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
    inline float cmf_y_pdf(float x) {
        float sx = x - 559.8692016601562f;
        float s = 23.981721878051758f;
        float a = -169.86920166015625f;
        float b = 270.13079833984375f;
        return logistic_pdf(sx, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
    }

    inline float cmf_y_sample(float u) {
        float s = 23.981721878051758f;
        float a = -169.86920166015625f;
        float b = 270.13079833984375f;
        float Pa = logistic_cdf(a, s);
        float Pb = logistic_cdf(b, s);
        return inverse_logistic_cdf(Pa + (Pb - Pa) * u, s) + 559.8692016601562f;
    }
}

namespace CIE_2015_10deg
{
    inline float asymmetric_gaussian(float x, float mean, float sigma, float a) {
        float denom = sigma + a * (x - mean);
        if (denom < 1.0e-15f) { denom = 1.0e-15f; }
        if (sigma * 2.0f < denom) { denom = sigma * 2.0f; }
        float k = (x - mean) / denom;
        return expf(-k * k);
    }
    inline float cmf_x(float x) {
        float a = 0.42f * asymmetric_gaussian(x, 445.5849609375f, 31.146467208862305f, 0.06435633450746536f);
        float b = 1.16f * asymmetric_gaussian(x, 594.5570068359375f, 48.602108001708984f, -0.04772702232003212f);
        return (a + b - a * b * 42.559776306152344f) * 1.0008530432426956f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.8383178710938f, 66.54190826416016f, -0.026492968201637268f) * 1.004315086937574f;
    }
    inline float cmf_z(float x) {
        return 2.146832f * asymmetric_gaussian(x, 445.9251708984375f, 30.91781997680664f, 0.08379141241312027f) * 0.9975112815948937f;
    }

    inline float logistic_pdf(float x, float s)
    {
        float k = expf(-fabsf(x) / s);
        return k / ((1.0 + k) * (1.0 + k) * s);
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
    inline float cmf_y_pdf(float x) {
        float sx = x - 554.270751953125f;
        float s = 26.879621505737305f;
        float a = -164.270751953125f;
        float b = 275.729248046875f;
        return logistic_pdf(sx, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
    }

    inline float cmf_y_sample(float u) {
        float s = 26.879621505737305f;
        float a = -164.270751953125f;
        float b = 275.729248046875f;
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

    //{
    //    double result = 0;
    //    int nSample = 1000000;
    //    for (int i = 0; i < nSample; i++)
    //    {
    //        float x = glm::mix(-5.0f, 5.0f, (float)i / nSample);
    //        float v = logistic_pdf(x, 1.0f);

    //        float dx = (float)10 / nSample;
    //        result += dx * v;
    //    }

    //    printf("");
    //}

    //double integral = 0;
    //for (int nm = 390; nm < 830; nm++)
    //{
    //    integral += cmf_y(nm);
    //}

    double integral_p = 0;
    for (int nm = 390; nm < 830; nm++)
    {
        integral_p += 1.0f / (830 - 390) * cmf_y_pdf(nm);
    }

    PCG rng;
    double sum = 0;
    int nSample = 10000;
    for (int i = 0; i < nSample; i++)
    {
        float lambda = cmf_y_sample(rng.uniformf());
        float p_lambda = cmf_y_pdf(lambda);
        float c = cmf_y(lambda) / p_lambda;
        sum += c;
    }
    double integral = sum / nSample;

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
        static float pdf_view_scale = 200;

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

        PrimBegin(PrimitiveMode::LineStrip, 2);
        for (int i = 0; i < N; i++)
        {
            float nm = glm::mix(360.0f, 830.0f, (float)i / N);
            float x = nm / 100.0f;
            float n = 1.64732f + 7907.16861f / (nm * nm);
            PrimVertex({ x, n * 10, 0 }, { 255, 255, 255 });
        }
        PrimEnd();

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
            ImGui::SliderFloat("pdf view scale", &pdf_view_scale, 0, 300);
        }
        ImGui::Checkbox("show Sampled Histogram(y)", &showSampledHistogram);

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
