import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.attitudes.FrameAlignedProvider;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.TideSystem;
import org.orekit.forces.gravity.potential.UnnormalizedSphericalHarmonicsProvider;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.EquinoctialOrbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngleType;
import org.orekit.propagation.PropagationType;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.semianalytical.dsst.DSSTPropagator;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTForceModel;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTJ2SquaredClosedForm;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTZonal;
import org.orekit.propagation.semianalytical.dsst.forces.ZeisModel;
import org.orekit.time.AbsoluteDate;

/** Offline reference: no UTC, EOP, gravity downloads, or celestial ephemerides. */
public final class LongHorizonReference {
    private static final double MU = 3.986004418e14;
    private static final double RADIUS = 6378137.0;
    private static final double J2 = 1.08262668e-3;

    private static List<DSSTForceModel> forces(String model) {
        UnnormalizedSphericalHarmonicsProvider gravity = GravityFieldFactory.getUnnormalizedProvider(
            RADIUS, MU, TideSystem.UNKNOWN, new double[][]{{1.0}, {0.0}, {-J2}},
            new double[][]{{0.0}, {0.0}, {0.0}});
        List<DSSTForceModel> result = new ArrayList<>();
        result.add(new DSSTZonal(gravity, 2, 1, 5));
        if (model.equals("j2_j2sq")) result.add(new DSSTJ2SquaredClosedForm(new ZeisModel(), gravity));
        return result;
    }

    private static void row(BufferedWriter writer, String model, String initialType,
                            double days, String outputType, SpacecraftState state) throws Exception {
        EquinoctialOrbit o = (EquinoctialOrbit) OrbitType.EQUINOCTIAL.convertType(state.getOrbit());
        Vector3D p = o.getPVCoordinates().getPosition();
        Vector3D v = o.getPVCoordinates().getVelocity();
        writer.write(String.format(Locale.ROOT,
            "%s,%s,%.1f,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g%n",
            model, initialType, days, outputType, o.getA(), o.getEquinoctialEx(), o.getEquinoctialEy(),
            o.getHx(), o.getHy(), o.getLM(), p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()));
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("Expected output CSV path");
        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(args[0]))) {
            writer.write("model,initial_type,days,output_type,a,ex,ey,hx,hy,lm,x,y,z,vx,vy,vz\n");
            for (String model : new String[]{"j2", "j2_j2sq"}) {
                for (String initialType : new String[]{"mean", "osculating"}) {
                    List<DSSTForceModel> models = forces(model);
                    FrameAlignedProvider attitude = new FrameAlignedProvider(FramesFactory.getGCRF());
                    EquinoctialOrbit initialOrbit = new EquinoctialOrbit(26560000.0, 0.01, -0.004, 0.4, 0.3, 0.8,
                        PositionAngleType.MEAN, FramesFactory.getGCRF(), AbsoluteDate.J2000_EPOCH, MU);
                    SpacecraftState initial = new SpacecraftState(initialOrbit);
                    if (initialType.equals("osculating")) {
                        initial = DSSTPropagator.computeMeanState(initial, attitude, models, 1e-14, 200);
                    }
                    // Reference tolerances tighter than the production C++
                    // default; one-hour cap independently bounds long-span error.
                    double[] absolute = {1e-6, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-12};
                    double[] relative = {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14};
                    DSSTPropagator propagator = new DSSTPropagator(
                        new DormandPrince853Integrator(1e-3, 3600.0, absolute, relative), PropagationType.MEAN, attitude);
                    propagator.setMu(MU);
                    for (DSSTForceModel force : models) propagator.addForceModel(force);
                    propagator.setInitialState(initial, PropagationType.MEAN);
                    for (double days : new double[]{0.0, 1.0, 30.0, 365.0}) {
                        SpacecraftState mean = days == 0.0 ? initial : propagator.propagate(AbsoluteDate.J2000_EPOCH.shiftedBy(days * 86400.0));
                        row(writer, model, initialType, days, "mean", mean);
                        SpacecraftState osculating = DSSTPropagator.computeOsculatingState(mean, attitude, forces(model));
                        row(writer, model, initialType, days, "osculating", osculating);
                    }
                }
            }
        }
    }
}
