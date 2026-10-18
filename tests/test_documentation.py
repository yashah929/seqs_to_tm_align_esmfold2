"""Validate requested local documentation/figure links and SVG source integrity."""
from pathlib import Path
import re
import unittest
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]


class DocumentationTests(unittest.TestCase):
    def test_relative_links_and_heading_anchors(self):
        files = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
        for source in files:
            for target in re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", source.read_text()):
                url = urlsplit(target)
                if url.scheme or url.netloc:
                    continue
                path = (source.parent / unquote(url.path)).resolve() if url.path else source
                with self.subTest(source=source.name, target=target):
                    self.assertTrue(path.exists(), f"Missing local link: {target}")
                    if url.fragment and path.suffix == ".md":
                        headings = re.findall(r"^#+\s+(.*)$", path.read_text(), re.M)
                        slugs = [re.sub(r"[^\w\- ]", "", h.lower()).replace(" ", "-") for h in headings]
                        self.assertIn(url.fragment, slugs)

    def test_figures_have_accessible_labels(self):
        namespace = {"s": "http://www.w3.org/2000/svg"}
        paths = sorted((ROOT / "docs/assets").glob("*.svg"))
        self.assertEqual(len(paths), 3)
        for path in paths:
            with self.subTest(path=path.name):
                svg = ET.parse(path).getroot()
                self.assertTrue(svg.find("s:title", namespace).text)
                self.assertTrue(svg.find("s:desc", namespace).text)
                self.assertEqual(svg.attrib["viewBox"].split()[2], "1000")


if __name__ == "__main__":
    unittest.main()
