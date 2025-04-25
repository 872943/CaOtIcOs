# Configura esto
$startIndex = 14
$folderPath = "C:\Users\USUARIO\Downloads\EJECUTABLES"
$total = 21

# Núcleos disponibles
$maxProcesses = [Environment]::ProcessorCount

# Genera el orden con condiciones periódicas de contorno
$order = for ($i = 0; $i -lt $total; $i++) {
    ($startIndex + $i) % $total
}

$running = @()
$executed = @()

for ($i = 0; $i -lt $order.Count; $i++) {
    $index = $order[$i]
    $exePath = Join-Path $folderPath "$index.exe"

    # Verifica si el archivo existe
    if (-Not (Test-Path $exePath)) {
        Write-Host "❌ No se encontró: $exePath" -ForegroundColor Red
        continue
    }

    # Espera si hay demasiados procesos activos
    while ($running.Count -ge $maxProcesses) {
        $running = $running | Where-Object { -not $_.HasExited }
        Start-Sleep -Milliseconds 200
    }

    # Muestra barra de progreso
    $percent = [math]::Round(($i / $total) * 100)
    Write-Progress -Activity "Ejecutando archivos .exe" -Status "$index.exe en progreso..." -PercentComplete $percent

    # Inicia proceso
    $proc = Start-Process -FilePath $exePath -PassThru
    $running += $proc
    $executed += "$index.exe"

    # Muestra lista actualizada de ejecutados
    Write-Host "`n✅ Ejecutados hasta ahora:" -ForegroundColor Green
    $executed | ForEach-Object { Write-Host "  $_" }
}

# Espera a que terminen todos
$running | ForEach-Object { $_.WaitForExit() }

Write-Host "`n🎉 Todos los programas que existían se han ejecutado." -ForegroundColor Cyan
