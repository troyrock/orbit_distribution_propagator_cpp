# Offline independent Java parity reference

From the simulator repository root, with an already built local Orekit 13.1.6
checkout and Hipparchus jars already in your Maven repository:

```powershell
.\tools\java_reference\generate.ps1 `
    -OrekitRoot D:\orekit\orekit-13.1.6 `
    -JavaHome 'C:\Program Files\Eclipse Adoptium\jdk-17.0.20.8-hotspot' `
    -MavenRepository "$env:USERPROFILE\.m2\repository"
```

`JavaHome` defaults to `JAVA_HOME`, or an installed JDK 17 if that variable
is unset. Use the current installation path on your machine. Maven itself
is not needed: the script reads the Hipparchus version from Orekit's POM and
uses existing jars. Missing prebuilt dependencies fail explicitly. The
script writes `tests/data/orekit_long_horizon.csv` and ignored generated
classes/libraries below `outputs/java-reference`.

The C++ probe accepts an output CSV and an optional reference CSV:

```powershell
.\build\distribution_backend_probe.exe `
    .\outputs\cpp_long_horizon.csv `
    .\tests\data\orekit_long_horizon.csv
```

The output parent directory must exist. A mismatch returns nonzero. Normal
CTest execution uses the committed fixture and requires no Java installation.
See [`tests/data/README.md`](../../tests/data/README.md) for force settings,
integration tolerances, hashes, measured discrepancies, and regression gates.
