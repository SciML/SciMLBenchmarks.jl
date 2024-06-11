within ;
package RC_Circuit

  model RC_Circuit_Base
    parameter Integer N=1000;
    Modelica.Electrical.Analog.Sources.SineVoltage source(f=0.05, V=1);
    Modelica.Electrical.Analog.Basic.Ground g;
    Modelica.Electrical.Analog.Basic.Resistor R[N](each R=1.0);
    Modelica.Electrical.Analog.Basic.Capacitor C[N](each C=1.0, each v(fixed=true));

  equation
    connect(g.p, source.n);
    connect(R[1].p, source.p);
    connect(C[1].p, R[1].n);
    connect(C[1].n, g.p);
    for i in 2:N loop
      connect(R[i].n, C[i].p);
      connect(R[i].p, R[i-1].n);
      connect(C[i].n, g.p);
    end for;
            annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end RC_Circuit_Base;

  package Test
    model RC_Circuit_MTK_test_5
      extends RC_Circuit_Base(N=5);
    end RC_Circuit_MTK_test_5;

    model RC_Circuit_MTK_test_10
      extends RC_Circuit_Base(N=10);
    end RC_Circuit_MTK_test_10;

    model RC_Circuit_MTK_test_20
      extends RC_Circuit_Base(N=20);
    end RC_Circuit_MTK_test_20;

    model RC_Circuit_MTK_test_40
      extends RC_Circuit_Base(N=40);
    end RC_Circuit_MTK_test_40;

    model RC_Circuit_MTK_test_60
      extends RC_Circuit_Base(N=60);
    end RC_Circuit_MTK_test_60;

    model RC_Circuit_MTK_test_80
      extends RC_Circuit_Base(N=80);
    end RC_Circuit_MTK_test_80;

    model RC_Circuit_MTK_test_160
      extends RC_Circuit_Base(N=160);
    end RC_Circuit_MTK_test_160;

    model RC_Circuit_MTK_test_320
      extends RC_Circuit_Base(N=320);
    end RC_Circuit_MTK_test_320;

    model RC_Circuit_MTK_test_480
      extends RC_Circuit_Base(N=480);
    end RC_Circuit_MTK_test_480;

    model RC_Circuit_MTK_test_640
      extends RC_Circuit_Base(N=640);
    end RC_Circuit_MTK_test_640;

    model RC_Circuit_MTK_test_800
      extends RC_Circuit_Base(N=800);
    end RC_Circuit_MTK_test_800;

    model RC_Circuit_MTK_test_1000
      extends RC_Circuit_Base(N=1000);
    end RC_Circuit_MTK_test_1000;

    model RC_Circuit_MTK_test_2000
      extends RC_Circuit_Base(N=2000);
    end RC_Circuit_MTK_test_2000;

    model RC_Circuit_MTK_test_3000
      extends RC_Circuit_Base(N=3000);
    end RC_Circuit_MTK_test_3000;

    model RC_Circuit_MTK_test_4000
      extends RC_Circuit_Base(N=4000);
    end RC_Circuit_MTK_test_4000;

    model RC_Circuit_MTK_test_5000
      extends RC_Circuit_Base(N=5000);
    end RC_Circuit_MTK_test_5000;

    model RC_Circuit_MTK_test_6000
      extends RC_Circuit_Base(N=6000);
    end RC_Circuit_MTK_test_6000;

    model RC_Circuit_MTK_test_7000
      extends RC_Circuit_Base(N=7000);
    end RC_Circuit_MTK_test_7000;

    model RC_Circuit_MTK_test_8000
      extends RC_Circuit_Base(N=8000);
    end RC_Circuit_MTK_test_8000;

    model RC_Circuit_MTK_test_9000
      extends RC_Circuit_Base(N=9000);
    end RC_Circuit_MTK_test_9000;

    model RC_Circuit_MTK_test_10000
      extends RC_Circuit_Base(N=10000);
    end RC_Circuit_MTK_test_10000;

    model RC_Circuit_MTK_test_20000
      extends RC_Circuit_Base(N=20000);
    end RC_Circuit_MTK_test_20000;
  end Test;

  package TimingFunctions

    function tic
      "Function to record the internal time in [s]"
      output Real tic;
    algorithm
      (ms_tic, sec_tic, min_tic, hour_tic, day_tic, mon_tic, year_tic) := Modelica.Utilities.System.getTime();
      tic := (ms_tic*0.001) + sec_tic + (min_tic*60) + (hour_tic*3600) + (day_tic*86400);
    end tic;

  end TimingFunctions;
    
end RC_Circuit;