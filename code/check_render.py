"""Optional: render site/index.html in headless Chromium and save a screenshot (needs `pip install playwright`)."""
import asyncio, sys
from pathlib import Path
from playwright.async_api import async_playwright

SITE = Path(__file__).resolve().parent.parent / "site" / "index.html"
OUT = Path(__file__).resolve().parent / "build" / "screenshot.png"

async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        pg = await b.new_page(viewport={"width": 1300, "height": 1000})
        errors = []
        pg.on("pageerror", lambda e: errors.append(str(e)))
        await pg.goto(SITE.as_uri())
        await pg.wait_for_timeout(8000)
        OUT.parent.mkdir(exist_ok=True)
        await pg.screenshot(path=str(OUT), full_page=True)
        await b.close()
    print(f"screenshot: {OUT}")
    if errors:
        print("page errors:", *errors, sep="\n  "); sys.exit(1)

asyncio.run(main())
