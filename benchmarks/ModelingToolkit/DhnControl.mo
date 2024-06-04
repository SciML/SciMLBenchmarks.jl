package DhnControl
  package Models
    connector FluidPort
      Real p;
      flow Real m;
      stream Real T;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {1, 132, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end FluidPort;

    connector HeatPort
      Real T;
      flow Real Q;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end HeatPort;

    model FluidRegion
      parameter Real T0 = 0.0;
      parameter Real lumped_T = 20;
      parameter Real L = 100;
      parameter Integer N = 100;
      parameter Real dn = 0.05;
      parameter Real eps = 1e-4;
      parameter Real C_shift = 0.0;
      parameter Real Rw = 0.0;
      parameter Real dx = L / N;
      parameter Real A = pi * dn ^ 2 / 4;
      //variables
      Real T[N](each start = T0);
      Real Twall[N](each start = T0);
      Real S[N];
      Real C[N];
      Real u;
      Real Re;
      Real Dxx;
      Real Pr;
      Real alpha;
      Real f;
      Real Nu;
      constant Real c[3] = {-1 / 8, -3 / 8, -3 / 8};
      constant Real pi = Modelica.Constants.pi;
      DhnControl.Models.FluidPort inlet annotation(
        Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.FluidPort outlet annotation(
        Placement(visible = true, transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.HeatPort heatport[N] annotation(
        Placement(visible = true, transformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      Re = 0.1 + dn * abs(u) / (1.7631408357540536e-6 - 5.37008081108929e-8 * lumped_T + 9.71612371740893e-10 * lumped_T ^ 2 - 1.0026133706457947e-11 * lumped_T ^ 3 + 5.3314727276417887e-14 * lumped_T ^ 4 - 1.1234839851121354e-16 * lumped_T ^ 5);
      Pr = 13.167521706890668 - 0.4441712774475796 * lumped_T + 0.008404050163010567 * lumped_T ^ 2 - 8.863470332920757e-5 * lumped_T ^ 3 + 4.768646349109186e-7 * lumped_T ^ 4 - 1.0118106644874892e-9 * lumped_T ^ 5;
      if Re <= 2300 then
        Nu = 3.66;
      elseif Re > 2300 and Re <= 3100 then
        Nu = 3.5239 * (Re / 1000) ^ 4 - 45.158 * (Re / 1000) ^ 3 + 212.13 * (Re / 1000) ^ 2 - 427.45 * (Re / 1000) + 316.08;
      else
        Nu = f / 8 * ((Re - 1000) * Pr) / (1 + 12.7 * (f / 8) ^ (1 / 2) * (Pr ^ (2 / 3) - 1));
      end if;
      if Re < 1000 then
        Dxx = dn ^ 2 / 4 * u ^ 2 / 48 / 0.14e-6;
      else
        Dxx = dn * u * (1.17e9 * Re ^ (-2.5) + 0.41);
      end if;
      f = 8 * ((8 / Re) ^ 12 + 1 / ((-2.457 * log((7 / Re) ^ 0.9 + 0.27 * (eps / dn))) ^ 16 + (37530 / Re) ^ 16) ^ 1.5) ^ (1 / 12);
//churchill
      alpha = Nu * (0.5629006705766969 + 0.002027906870916878 * lumped_T - 1.0062563416770859e-5 * lumped_T ^ 2 + 1.2897392253800086e-8 * lumped_T ^ 3) / dn;
      inlet.m = -outlet.m;
      inlet.p = outlet.p;
      inlet.T = inStream(inlet.T);
      outlet.T = T[N];
      u = inlet.m / (1000.2325951116342 + 0.02341353916404164 * inlet.T - 0.006409744347791998 * inlet.T ^ 2 + 2.6080324835019334e-5 * inlet.T ^ 3 - 6.044425395404277e-8 * inlet.T ^ 4) / A;
      C[:] = {dx * A * (4.2146276665115245e6 - 2231.4603937980637 * T[i] + 24.383955576707972 * T[i] ^ 2 - 0.3896654168053555 * T[i] ^ 3 + 0.002537945351479517 * T[i] ^ 4 - 5.879294024916786e-6 * T[i] ^ 5) for i in 1:N};
      S[:] = {heatport[i].Q for i in 1:N};
      Twall[:] = {heatport[i].T for i in 1:N};
      S[:] = {1 / (1 / (alpha * dn * pi * dx) + abs(Rw / 1000)) * (Twall[i] - T[i]) for i in 1:N};
      der(T[1]) = u / dx * (inlet.T - T[1]) + Dxx * (T[2] - T[1]) / dx ^ 2 + S[1] / (C[1] - C_shift);
      der(T[2]) = u / dx * (c[1] * inlet.T - sum(c) * T[1] + c[2] * T[2] + c[3] * T[3]) + Dxx * (T[1] - 2 * T[2] + T[3]) / dx ^ 2 + S[2] / (C[2] - C_shift);
      der(T[3:N - 1]) = {u / dx * (c[1] * T[i - 2] - sum(c) * T[i - 1] + c[2] * T[i] + c[3] * T[i + 1]) + Dxx * (T[i - 1] - 2 * T[i] + T[i + 1]) / dx ^ 2 + S[i] / (C[i] - C_shift) for i in 3:N - 1};
      der(T[N]) = u / dx * (T[N - 1] - T[N]) + Dxx * (T[N - 1] - T[N]) / dx ^ 2 + S[N] / (C[N] - C_shift);
      annotation(
        Icon(graphics = {Rectangle(fillColor = {50, 151, 213}, fillPattern = FillPattern.Solid, extent = {{-80, 20}, {80, -20}})}));
    end FluidRegion;

    model Source
      parameter Real t_start = 12 * 3600;
      parameter Real p_feed = 100000;
      parameter Real m_flow = 1.0;
      parameter Real T_out = 1.0;
      FluidPort outlet annotation(
        Placement(visible = true, transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      if time > t_start then
        outlet.T = 56.0 + 12.0;
      else
        outlet.T = 12.0;
      end if;
      outlet.m = -m_flow;
      outlet.p = p_feed;
      annotation(
        Icon(graphics = {Rectangle(origin = {-10, 0}, fillColor = {204, 204, 204}, fillPattern = FillPattern.Solid, extent = {{-90, 100}, {90, -100}}), Text(origin = {-12, 52}, extent = {{-70, 24}, {70, -24}}, textString = "%name")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end Source;

    model Sink
      FluidPort inlet annotation(
        Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      inlet.T = inStream(inlet.T);
      annotation(
        Icon(graphics = {Rectangle(origin = {10, 0}, fillColor = {212, 212, 212}, fillPattern = FillPattern.Solid, extent = {{-90, 100}, {90, -100}}), Text(origin = {3, 53}, extent = {{-63, 25}, {63, -25}}, textString = "%name")}));
    end Sink;

    model CylindricalWallFem
      parameter Real T0 = 0.0;
      parameter Real d = 0.3;
      parameter Real L = 100;
      parameter Integer N = 100;
      parameter Real dx = L / N;
      parameter Real t_layer[2] = {0.00391, 0.013};
      parameter Integer Ne = size(t_layer,1);
      parameter Real cp[Ne] = {500, 1200};
      parameter Real rho[Ne] = {7850, 40};
      parameter Real lambda[Ne] = {50, 0.04};
      parameter Integer Nn = Ne + 1;
      parameter Real dn[Nn] = cat(1, {d}, d .+ 2.0.*{sum(t_layer[1:i]) for i in 1:Ne});
      parameter Real Cn[Nn] = cat(1, {pi/8*(dn[i+1]^2-dn[i]^2)*cp[i]*rho[i] for i in 1:Ne}, {0})+ cat(1, {0}, {pi/8*(dn[i+1]^2-dn[i]^2)*cp[i]*rho[i] for i in 1:Ne});
      parameter Real Ke[Ne] = {2*pi*lambda[i]/log(dn[i+1]/dn[i]) for i in 1:Ne};
      Real Tn[N,Nn](each start = T0);
      Real Qe[N,Ne];
      DhnControl.Models.HeatPort inner_heatport[N] annotation(
          Placement(visible = true, transformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.HeatPort outer_heatport[N] annotation(
          Placement(visible = true, transformation(origin = {0, 92}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      constant Real pi = Modelica.Constants.pi;
    equation
      for i in 1:N loop
        inner_heatport[i].T = Tn[i,1];
        outer_heatport[i].T = Tn[i,Nn];
        Qe[i,:] = {Ke[j]*(-Tn[i,j+1]+Tn[i,j]) for j in 1:Ne};
      end for;

      der(Tn[:,1])*(Cn[1]) = {inner_heatport[i].Q - Qe[i,1] for i in 1:N};
      for i in 1:N loop
        der(Tn[i,2:Ne]).*Cn[2:Ne] = {Qe[i,j-1]-Qe[i,j] for j in 2:Ne};
      end for;
      der(Tn[:,Nn])*Cn[Nn] = {Qe[i, Ne] + outer_heatport[i].Q for i in 1:N};

      annotation(
        Icon(graphics = {Rectangle(fillColor = {214, 197, 109}, fillPattern = FillPattern.Solid, extent = {{-100, 80}, {100, -80}}), Text(origin = {-4, 9}, extent = {{-72, 27}, {72, -27}}, textString = "%name")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end CylindricalWallFem;

    model PreinsulatedPipe
      parameter Real T0 = 0.0;
      parameter Real lumped_T = 20.0;
      parameter Real eps = 1e-4;
      parameter Real L =470;
      parameter Integer N = 100;
      parameter Real dn = 0.3;
      parameter Real t_layer[:] = {0.0056, 0.058};
      parameter Real cp[:] = {500, 1200};
      parameter Real rho[:] = {7850, 40};
      parameter Real lambda[:] = {50, 0.04};
      parameter Real alpha_out = 4.0;
      parameter Real Tenv = 18.0;
      DhnControl.Models.FluidPort inlet annotation(
        Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.FluidPort outlet annotation(
        Placement(visible = true, transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.FluidRegion fluidRegion(L=L, N=N, dn=dn, T0=T0, lumped_T=lumped_T,eps=eps) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.CylindricalWallFem cylindricalWallFem(t_layer=t_layer, L=L, N=N, d=dn, rho=rho, lambda=lambda, cp=cp, T0=T0) annotation(
        Placement(visible = true, transformation(origin = {0, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DhnControl.Models.CylindricalSurfaceConvection cylindricalSurfaceConvection(L=L, N=N, d=dn+2*sum(t_layer), alpha=alpha_out,Tenv=Tenv) annotation(
        Placement(visible = true, transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(fluidRegion.outlet, outlet) annotation(
        Line(points = {{9, 0}, {90, 0}}));
      connect(fluidRegion.inlet, inlet) annotation(
        Line(points = {{-9, 0}, {-90, 0}}));
      connect(fluidRegion.heatport, cylindricalWallFem.inner_heatport) annotation(
        Line(points = {{0, 4}, {0, 16}}, thickness = 0.5));
  connect(cylindricalWallFem.outer_heatport, cylindricalSurfaceConvection.heatport) annotation(
        Line(points = {{0, 34}, {0, 42}}, thickness = 0.5));
    annotation(
        Icon(graphics = {Rectangle( lineColor = {254, 247, 33}, fillColor = {108, 115, 77}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-80, 34}, {80, -34}}),Rectangle(fillColor = {108, 108, 108}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-80, 24}, {80, -24}}), Rectangle(lineColor = {206, 206, 206},fillColor = {26, 139, 238}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-80, 20}, {80, -20}})}));
        end PreinsulatedPipe;

    model CylindricalSurfaceConvection
    parameter Real d = 0.3;
    parameter Real L = 470;
    parameter Integer N = 100;
    parameter Real Tenv = 18;
    parameter Real dx= L/N;
    parameter Real alpha = 4.0;
    constant Real pi = Modelica.Constants.pi;
    DhnControl.Models.HeatPort heatport[N] annotation(
        Placement(visible = true, transformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
for i in 1:N loop
    heatport[i].Q = alpha*pi*d*dx*(heatport[i].T-Tenv);
    end for;
    annotation(
        Icon(graphics = {Line(origin = {-1.45, -60.62}, points = {{-100, 0}, {100, 0}}, thickness = 1), Line(origin = {-70.06, 0.679544}, points = {{9.76713, -58.9684}, {-10.2329, -20.9684}, {9.76713, 11.0316}, {-8.23287, 59.0316}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Line(origin = {-30.7239, 0.513569}, points = {{9.76713, -58.9684}, {-10.2329, -20.9684}, {9.76713, 11.0316}, {-8.23287, 59.0316}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Line(origin = {10.9359, 0.63805}, points = {{9.76713, -58.9684}, {-10.2329, -20.9684}, {9.76713, 11.0316}, {-8.23287, 59.0316}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Line(origin = {52.5957, 0.762531}, points = {{9.76713, -58.9684}, {-10.2329, -20.9684}, {9.76713, 11.0316}, {-8.23287, 59.0316}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier)}));
        end CylindricalSurfaceConvection;
  end Models;

  package Test
    model test_adiabatic_470
      DhnControl.Models.Source source(m_flow = 2.75) annotation(
        Placement(visible = true, transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DhnControl.Models.FluidRegion fluidRegion(L = 470, N = 470, T0 = 12, dn = 0.3127) annotation(
        Placement(visible = true, transformation(origin = {-12, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Models.Sink sink annotation(
        Placement(visible = true, transformation(origin = {28, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, fluidRegion.inlet) annotation(
        Line(points = {{-40, 10}, {-21, 10}}));
      connect(fluidRegion.outlet, sink.inlet) annotation(
        Line(points = {{-3, 10}, {20, 10}}));
      annotation(
        experiment(StartTime = 0, StopTime = 68400, Tolerance = 1e-06, Interval = 100));
    end test_adiabatic_470;

    model test_preinsulated_470
      DhnControl.Models.Source source(m_flow = 2.75) annotation(
        Placement(visible = true, transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Models.Sink sink annotation(
        Placement(visible = true, transformation(origin = {28, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Models.PreinsulatedPipe preinsulatedPipe(L = 470, N = 1280, T0 = 12, dn = 0.3127)  annotation(
        Placement(visible = true, transformation(origin = {-12, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(source.outlet, preinsulatedPipe.inlet) annotation(
        Line(points = {{-40, 10}, {-20, 10}}));
    connect(preinsulatedPipe.outlet, sink.inlet) annotation(
        Line(points = {{-2, 10}, {20, 10}}));
      annotation(
        experiment(StartTime = 0, StopTime = 68400, Tolerance = 1e-06, Interval = 100),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "ida"));
    end test_preinsulated_470;

    model test_preinsulated_470_5
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=5));
    end test_preinsulated_470_5;

    model test_preinsulated_470_10
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=10));
    end test_preinsulated_470_10;

    model test_preinsulated_470_20
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=20));
    end test_preinsulated_470_20;

    model test_preinsulated_470_40
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=40));
    end test_preinsulated_470_40;

    model test_preinsulated_470_60
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=60));
    end test_preinsulated_470_60;

    model test_preinsulated_470_80
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=80));
    end test_preinsulated_470_80;

    model test_preinsulated_470_160
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=160));
    end test_preinsulated_470_160;

    model test_preinsulated_470_320
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=320));
    end test_preinsulated_470_320;

    model test_preinsulated_470_480
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=480));
    end test_preinsulated_470_480;

    model test_preinsulated_470_640
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=640));
    end test_preinsulated_470_640;

    model test_preinsulated_470_800
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=800));
    end test_preinsulated_470_800;

    model test_preinsulated_470_960
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=960));
    end test_preinsulated_470_960;

    model test_preinsulated_470_1280
        extends DhnControl.Test.test_preinsulated_470(preinsulatedPipe(N=1280));
    end test_preinsulated_470_1280;

  end Test;
end DhnControl;
