## A simple and table-free analytic approximation of color matching functions with importance sampling

- CIE 1931 2 Degree 
- CIE 2015 10 Degree Standard Observer

The curves have slightly asymmetric shape. Thus, I used a modified version gaussian function:

![image](AG.gif)

And fitted by PyTorch.

Because it's hard to find a simple CDF, I used [Trimmed Logistic](https://pbr-book.org/4ed/Utilities/Mathematical_Infrastructure#TrimmedLogistic) for importance sampling.

## CIE 2015 10 Degree Standard Observer
![image](2015.png)

```
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
        return a + b - a * b * 42.559776306152344f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.8383178710938f, 66.54190826416016f, -0.026492968201637268f);
    }
    inline float cmf_z(float x) {
        return 2.146832f * asymmetric_gaussian(x, 445.9251708984375f, 30.91781997680664f, 0.08379141241312027f);
    }
}

```

<details> <summary>Importance Sampling( cmf_y_pdf( nm ), cmf_y_sample( 0.0 to 1.0 ) )</summary>

```
namespace CIE_2015_10deg
{
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
```
</details> 

## CIE 1931 2 Degree 
![image](1931.png)

```
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
        return a + b - a * b * 21.016616821289062f;
    }
    inline float cmf_y(float x) {
        return 1.0f * asymmetric_gaussian(x, 556.5616455078125f, 59.5950927734375f, 0.056370146572589874f);
    }
    inline float cmf_z(float x) {
        return 1.7829682f * asymmetric_gaussian(x, 447.90643310546875f, 32.452659606933594f, 0.12635648250579834f);
    }
}

```

<details> <summary>Importance Sampling( cmf_y_pdf( nm ), cmf_y_sample( 0.0 to 1.0 ) )</summary>

```
namespace CIE_1931_2deg
{
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
```

</details>

## Build

```
git submodule update --init
premake5 vs2022
```

## Ref
- [Simple Analytic Approximations to the CIE XYZ Color Matching Functions](https://jcgt.org/published/0002/02/01/)
- [colour-science](https://colour.readthedocs.io/en/v0.4.5/generated/colour.MSDS_CMFS.html)
- [PBRT Trimmed Logistic](https://pbr-book.org/4ed/Utilities/Mathematical_Infrastructure#TrimmedLogistic)

