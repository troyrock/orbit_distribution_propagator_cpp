param(
    [string]$OrekitRoot = 'D:\orekit\orekit-13.1.6',
    [string]$JavaHome = $env:JAVA_HOME,
    [string]$MavenRepository = (Join-Path $env:USERPROFILE '.m2\repository'),
    [string]$Output = (Join-Path $PSScriptRoot '..\..\tests\data\orekit_long_horizon.csv')
)
$ErrorActionPreference = 'Stop'
if (-not $JavaHome) {
    $JavaHome = (Get-ChildItem -LiteralPath (Join-Path $env:ProgramFiles 'Eclipse Adoptium') -Directory |
        Where-Object { $_.Name -like 'jdk-17*' } | Sort-Object Name -Descending | Select-Object -First 1).FullName
}
$java = Join-Path $JavaHome 'bin\java.exe'
$javac = Join-Path $JavaHome 'bin\javac.exe'
foreach ($tool in @($java, $javac)) {
    if (-not (Test-Path -LiteralPath $tool)) { throw "Required JDK tool missing: $tool" }
}
[xml]$pom = Get-Content -LiteralPath (Join-Path $OrekitRoot 'pom.xml') -Raw
$version = $pom.project.properties.'orekit.hipparchus.version'
$entries = @((Join-Path $OrekitRoot 'target\orekit-13.1.6.jar'))
foreach ($component in @('core', 'geometry', 'ode', 'fitting', 'optim', 'filtering', 'stat')) {
    $entries += Join-Path $MavenRepository "org\hipparchus\hipparchus-$component\$version\hipparchus-$component-$version.jar"
}
foreach ($entry in $entries) {
    if (-not (Test-Path -LiteralPath $entry)) { throw "Required prebuilt offline dependency missing: $entry" }
}
$build = Join-Path $PSScriptRoot '..\..\outputs\java-reference\classes'
New-Item -ItemType Directory -Force -Path $build | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Output) | Out-Null
$localLibraries = Join-Path $PSScriptRoot '..\..\outputs\java-reference\libraries'
New-Item -ItemType Directory -Force -Path $localLibraries | Out-Null
$localEntries = @()
foreach ($entry in $entries) {
    # Some Windows restricted accounts can read Maven jars but Java ZIP close
    # cannot resolve their original real path. Workspace copies avoid that
    # runtime issue while keeping the upstream checkout and Maven cache intact.
    $localJar = Join-Path $localLibraries (Split-Path -Leaf $entry)
    Copy-Item -LiteralPath $entry -Destination $localJar -Force
    $localEntries += $localJar
}
$classpath = $localEntries -join [IO.Path]::PathSeparator
& $javac -cp $classpath -d $build (Join-Path $PSScriptRoot 'LongHorizonReference.java')
if ($LASTEXITCODE -ne 0) { throw 'Java reference compilation failed' }
& $java -cp ($build + [IO.Path]::PathSeparator + $classpath) LongHorizonReference $Output
if ($LASTEXITCODE -ne 0) { throw 'Java reference generation failed' }
Write-Output "Wrote $Output"
