import torch
import torch.nn as nn
import csv
import matplotlib.pyplot as plt
import colour # pip install colour-science


CMF_Name = 'cie_2_1931'
X_PEAK_1 = 0.42
X_PEAK_2 = 1.64

# CMF_Name = 'CIE 2015 10 Degree Standard Observer'
# X_PEAK_1 = 0.42
# X_PEAK_2 = 1.160

CMF: colour.colorimetry.XYZ_ColourMatchingFunctions = colour.MSDS_CMFS[CMF_Name]

wavelength = CMF.wavelengths
Xs = [xyz[0] for xyz in CMF.values]
Ys = [xyz[1] for xyz in CMF.values]
Zs = [xyz[2] for xyz in CMF.values]

# with open( f"{CMF_Name}.csv", 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile, delimiter=',')
#     for i in range(len(wavelength)):
#         writer.writerow([wavelength[i], Xs[i], Ys[i], Zs[i]])

# wavelength = []
# Xs = []
# Ys = []
# Zs = []
# with open(f"{CMF_Name}.csv", newline="", encoding="utf-8") as f:
#     reader = csv.reader(f)
#     for row in reader:
#         wavelength.append(float(row[0]))
#         Xs.append(float(row[1]))
#         Ys.append(float(row[2]))
#         Zs.append(float(row[3]))

def asymmetric_gaussian( x, mean, sigma, a ):
    k = (x - mean) / ( sigma + a * (x - mean) )
    return torch.exp( - k * k )

class AGaussianSinglePeak(nn.Module):
    def __init__(self, mean_init, peak):
        super().__init__()
        self.peak = peak
        self.mean = nn.Parameter(torch.tensor(mean_init))
        self.sigma = nn.Parameter(torch.tensor(100.0))
        self.a = nn.Parameter(torch.tensor(0.0))

    def forward(self, x):
        return self.peak * asymmetric_gaussian(x, self.mean, self.sigma, self.a)
    
    def cFunction(self, name):
        return f"""
        inline float asymmetric_gaussian( float x, float mean, float sigma, float a ) {{
            float k = (x - mean) / ( sigma + a * (x - mean) );
            return expf( - k * k );
        }}
        inline float {name}( float x ) {{
            return {self.peak}f * asymmetric_gaussian( x, {self.mean.item()}f, {self.sigma.item()}f, {self.a.item()}f );
        }}
        """ 
    
class AGaussianDualPeak(nn.Module):
    def __init__(self):
        super().__init__()
        self.f1 = AGaussianSinglePeak(mean_init=430.0, peak=X_PEAK_1)  
        self.f2 = AGaussianSinglePeak(mean_init=600.0, peak=X_PEAK_2)
        self.c = nn.Parameter(torch.tensor(10.0))
    def forward(self, x):
        a = self.f1(x)
        b = self.f2(x)
        return a + b - a * b * self.c
    
    def cFunction(self, name):
        return f"""
        inline float asymmetric_gaussian( float x, float mean, float sigma, float a ) {{
            float k = (x - mean) / ( sigma + a * (x - mean) );
            return expf( - k * k );
        }}
        inline float {name}( float x ) {{
            float a = {self.f1.peak}f * asymmetric_gaussian( x, {self.f1.mean.item()}f, {self.f1.sigma.item()}f, {self.f1.a.item()}f );
            float b = {self.f2.peak}f * asymmetric_gaussian( x, {self.f2.mean.item()}f, {self.f2.sigma.item()}f, {self.f2.a.item()}f );
            return a + b - a * b * {self.c.item()}f;
        }}
        """

def logistic( x, s ):
    k = torch.exp( - torch.abs(x) / s )
    return s * k / ( (1.0 + k) ** 2 )

class LogisticPDFSinglePeak(nn.Module):
    def __init__(self, mean_init, lowerBound, upperBound ):
        super().__init__()
        self.sigma = nn.Parameter(torch.tensor(16.0))
        self.mean = nn.Parameter(torch.tensor(mean_init))
        self.lowerBound = lowerBound
        self.upperBound = upperBound

    def forward(self, x):
        return logistic( x - self.mean, self.sigma )
    
    def cFunction(self, name_pdf, name_sample):
        return f"""
        inline float logistic_pdf( float x, float s )
        {{
            float k = expf( - fabsf(x) / s );
            return s * k / ( (1.0 + k) * (1.0 + k) );
        }}
        inline float logistic_cdf( float x, float s )
        {{
            return 1.0f / ( 1.0f + expf( - x / s ) );
        }}
        inline float inverse_logistic_cdf( float u, float s ) 
        {{
            if (0.99999994f < u) {{ u = 0.99999994f; }}
            if (u < 1.175494351e-38f) {{ u = 1.175494351e-38f; }}
            return -s * logf( 1.0f / u - 1.0f );
        }}
        inline float trimmed_logistic_pdf( float x, float s, float a, float b )
        {{
            return logistic_pdf( x, s ) / ( logistic_cdf(b, s) - logistic_cdf(a, s) );
        }}

        inline float {name_pdf}( float x ) {{
            float sx = x - {self.mean.item()}f;
            float s = {self.sigma.item()}f;
            float a = {self.lowerBound - self.mean.item()};
            float b = {self.upperBound - self.mean.item()};
            return logistic_pdf(sx, {self.sigma.item()}f) / ( logistic_cdf(b, s) - logistic_cdf(a, s) );
        }}

        inline float {name_sample}( float u ) {{
            float s = {self.sigma.item()}f;
            float a = {self.lowerBound - self.mean.item()};
            float b = {self.upperBound - self.mean.item()};
            float Pa = logistic_cdf(a, s);
            float Pb = logistic_cdf(b, s);
            return inverse_logistic_cdf( Pa + (Pb - Pa) * u, s ) + {self.mean.item()}f;
        }}
        """ 

if True:
    function_x = AGaussianDualPeak()
    optimizer = torch.optim.Adam(params=function_x.parameters(), lr=0.1, eps=1e-15)
    loss_fn = torch.nn.MSELoss()

    for i in range(20000):
        pred = function_x(torch.tensor(wavelength))

        loss = loss_fn(pred, torch.tensor(Xs))
        
        optimizer.zero_grad()
        
        loss.backward()

        optimizer.step()

    print(function_x.cFunction("cmf_x"))

    plt.plot(wavelength, Xs, label='X(reference)')
    plt.plot(wavelength, pred.detach().numpy(), label='X(fit)')
    plt.xlabel('Wavelength (nm)')
    plt.legend()
    plt.show()

if True:
    function_y = AGaussianSinglePeak(mean_init=550.0, peak=1.0)
    optimizer = torch.optim.Adam(params=function_y.parameters(), lr=0.1, eps=1e-15)
    loss_fn = torch.nn.MSELoss()

    for i in range(10000):
        pred1 = function_y(torch.tensor(wavelength))

        loss = loss_fn(pred1, torch.tensor(Ys))
        
        optimizer.zero_grad()
        
        loss.backward()

        optimizer.step()

    print(function_y.cFunction("cmf_y"))

    plt.plot(wavelength, Ys, label='Y(reference)')
    plt.plot(wavelength, pred1.detach().numpy(), label='Y(fit)')
    plt.xlabel('Wavelength (nm)')
    plt.legend()
    plt.show()

if True:
    function_z = AGaussianSinglePeak(mean_init=450.0, peak=max(Zs))
    optimizer = torch.optim.Adam(params=function_z.parameters(), lr=0.1, eps=1e-15)
    loss_fn = torch.nn.MSELoss()

    for i in range(10000):
        pred = function_z(torch.tensor(wavelength))

        loss = loss_fn(pred, torch.tensor(Zs))
        
        optimizer.zero_grad()
        
        loss.backward()

        optimizer.step()

    print(function_z.cFunction("cmf_z"))

    plt.plot(wavelength, Zs, label='Z(reference)')
    plt.plot(wavelength, pred.detach().numpy(), label='Z(fit)')
    plt.xlabel('Wavelength (nm)')
    plt.legend()
    plt.show()

if True:
    function_y_pdf = LogisticPDFSinglePeak(mean_init=550., lowerBound=390, upperBound=830)
    optimizer = torch.optim.Adam(params=function_y_pdf.parameters(), lr=0.1, eps=1e-15)
    loss_fn = torch.nn.KLDivLoss()

    for i in range(10000):
        pred = function_y_pdf(torch.tensor(wavelength))

        loss = loss_fn((pred / pred.sum()).log(), torch.tensor(Ys) / sum(Ys))
        
        optimizer.zero_grad()
        
        loss.backward()

        optimizer.step()

    print(function_y_pdf.cFunction("cmf_y_pdf", "cmf_y_sample"))

    peak = function_y_pdf(function_y_pdf.mean)
    plt.plot(wavelength, Ys, label='Y(reference)')
    plt.plot(wavelength, pred.detach().numpy() / peak.item(), label='Y(fit)')
    plt.xlabel('Wavelength (nm)')
    plt.legend()
    plt.show()

