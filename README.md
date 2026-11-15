# Kibon 3D

Modelo 3D interativo do Kibon a partir do levantamento com drone: MDT, ortofoto e curvas de nível de 1 m.

**Acesse:** https://ghagui.github.io/Kibon-3D/

## O que dá para ver

- **Terreno** com ortofoto ou cores por altitude, curvas de nível e exagero vertical.
- **Clima:** céu com sol que segue a hora do dia, nuvens, chuva com vento e relâmpagos.
- **Água:** rede de drenagem calculada sobre o MDT. Com chuva o solo encharca, os canais enchem aos poucos, as poças aparecem nas depressões e partículas mostram a água escorrendo até as saídas.
- **Risco de deslizamento:** mapa de suscetibilidade (chuva crítica de cada ponto, modelo de encosta infinita SHALSTAB) e simulação: com a chuva o solo encharca, as áreas instáveis piscam, as encostas cedem e a lama desce pelos caminhos da água até parar. Profundidade do solo, coesão e atrito são ajustáveis. É um modelo de estudo com parâmetros genéricos e não substitui o mapeamento oficial de risco (Defesa Civil / IPT).
- **Áreas verdes para retenção:** manchas de vegetação plana que recebem enxurrada de telhados e asfalto, com o volume que cada uma poderia segurar numa chuva de projeto (sem contar duas vezes a água já retida rio acima). É um potencial teórico que ignora bueiros e galerias.
- **Noite:** iluminação pública ao longo das ruas do OpenStreetMap e luzes nas casas, que acendem e apagam conforme a hora.
- **Inspeção:** passe o mouse (ou toque, no celular) para ver altitude, coordenadas, declividade, área de contribuição e fator de segurança da encosta.

No celular o visualizador usa automaticamente um perfil leve (malha e textura menores, menos partículas). Dá para trocar em *Terreno → Qualidade gráfica* ou com `?q=leve` / `?q=alta` na URL.

## Estrutura

| Pasta | Conteúdo |
|---|---|
| `index.html` | Visualizador (three.js, sem build) |
| `assets/` | Dados prontos para a web: MDT 1024 px com drenagem, curvas, postes, luzes e ortofoto 4k/2k |
| `dados/` | Dados geográficos otimizados: MDT 0,5 m e ortofoto 0,2 m (Cloud Optimized GeoTIFF), curvas em shapefile, vias e edificações do OSM |
| `tools/` | Scripts que geram `dados/` e `assets/` |

### Dados

- `MDT_kibon.tif`: 0,5 m, float32, DEFLATE (original de ~950 MB com 7 cm).
- `Ortofoto_Kibon.tif`: 0,2 m, JPEG com máscara (original de ~4,5 GB com 3 cm).
- `Curvas_1m.shp/.shx/.dbf/.prj`: SIRGAS 2000 / UTM 23S (EPSG:31983). O `.prj` e o `.shx` não vinham no original; o sistema foi confirmado pela sobreposição com o MDT.
- Fora da área do levantamento o relevo é estimado e o chão é gerado, só para dar contexto visual.
- Luzes das casas: o OSM tem pouquíssimas edificações mapeadas aqui, então os telhados foram detectados na ortofoto.

## Regerar

```bash
pip install -r tools/requirements.txt
# a partir dos originais em alta resolução
python tools/otimizar_dados.py <pasta_originais> dados
# assets do visualizador
python tools/build_assets.py dados assets
```

Para ver localmente: `python -m http.server` e abra http://localhost:8000.

O deploy no GitHub Pages roda sozinho a cada push na `main` (`.github/workflows/pages.yml`).

## Créditos

Vias e edificações © colaboradores do [OpenStreetMap](https://www.openstreetmap.org/copyright), licença ODbL.
