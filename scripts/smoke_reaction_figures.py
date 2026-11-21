"""Exercise scheme orientation, exports, inspector sizing and precursor changes through the UI."""
import argparse
import time
from pathlib import Path
from xml.etree import ElementTree as ET

from playwright.sync_api import sync_playwright


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-url', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with sync_playwright() as pw:
        chrome=Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe')
        browser=pw.chromium.launch(executable_path=str(chrome) if chrome.exists() else None)
        page=browser.new_page(viewport={'width':1450,'height':970},accept_downloads=True)
        page.set_default_timeout(60000)
        errors=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        def wait(js):
            until=time.monotonic()+60
            while not page.evaluate(js):
                if time.monotonic()>until: raise AssertionError((js,page.locator('#notice').text_content(),errors))
                page.wait_for_timeout(100)
        def save(name,feeds,ext='svg',crystal=False):
            wait("document.querySelector('.showing-scheme .reaction-figure')?.naturalWidth>0")
            button=page.locator('#save-crystal-drawing') if crystal else page.locator('#exports').get_by_role('button',name=ext.upper(),exact=True)
            with page.expect_download() as out: button.click()
            target=args.output/(name+'.'+ext)
            out.value.save_as(target)
            if ext=='svg':
                svg=ET.parse(target).getroot()
                assert len([k for k in svg if k.get('data-part','').startswith('precursor-')])==feeds
                assert len(svg.findall("./*[@data-part='framework']"))==1
                return float(svg.get('width')),float(svg.get('height'))
            assert target.read_bytes()[:8]==b'\x89PNG\r\n\x1a\n'
        page.goto(args.base_url)
        wait("window.COFdrawStudio?.state.online && document.querySelector('#build:not(:disabled)')")
        page.locator('#build').click()
        wait("window.COFdrawExtent?.visible && !document.querySelector('#save-project').disabled")
        assert page.locator('#build-periodic-from-pair').is_hidden()
        assert not page.locator('#source-select').count()
        page.locator('#appearance-options summary').first.click()
        for width in (1450,1050):
            page.set_viewport_size({'width':width,'height':970})
            boxes=page.locator('#style-list').evaluate("""s=>[s.clientWidth,s.scrollWidth,...[...s.children].map(b=>b.getBoundingClientRect().height)]""")
            assert boxes[0]>150 and boxes[1]<=boxes[0]+1 and max(boxes[2:])<55,boxes
        page.set_viewport_size({'width':1450,'height':970})
        page.locator('#appearance-options summary').first.click()
        btn=page.locator('.canvas-panel .scheme-toggle')
        btn.click()
        w,h=save('binary-horizontal',2)
        assert w>h
        page.locator('.canvas-panel .scheme-direction').select_option('vertical')
        wait("document.querySelector('.reaction-figure')?.naturalHeight>document.querySelector('.reaction-figure')?.naturalWidth")
        w,h=save('binary-vertical',2);assert h>w
        save('binary-vertical',2,'png')
        page.locator('#extent-options button[value="three-pores"]').click()
        wait("COFdrawExtent.visible?.summary?.pore_count===3")
        save('three-connected-pores',2)
        page.select_option('#assembly-count','1')
        wait("COFdrawExtent.visible&&COFdrawStudio.state.recipe?.assembly?.monomer_ids?.length===1")
        save('one-precursor',1)
        page.select_option('#assembly-count','3')
        wait("COFdrawExtent.visible&&COFdrawStudio.state.recipe?.assembly?.monomer_ids?.length===3")
        save('three-precursors',3)
        page.select_option('#assembly-count','2')
        wait("!document.querySelector('#cof-preset').disabled")
        page.select_option('#cof-preset','crystal:cof-300-c5')
        wait("window.COFdrawCrystalDesign?.report && document.querySelector('#crystal-drawing-image')?.naturalWidth>0")
        save('diamond-scheme',2,crystal=True)
        page.select_option('#crystal-drawing-view','Depth fade · two building blocks')
        save('diamond-depth',2,crystal=True)
        page.locator('#rotate-crystal-design').click()
        wait("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        page.locator('#close-science').click()
        assert page.locator('.crystal-drawing-page.showing-scheme .reaction-figure').is_visible()
        page.locator('#leave-crystal-design').click()
        wait("!COFdrawCrystalDesign.active && !document.querySelector('#build-periodic-from-pair').hidden")
        page.locator('#build-periodic-from-pair').click()
        wait("COFdrawCrystalDesign.report?.structure?.provenance?.network_dimension===3")
        save('custom-diamond',2,crystal=True)
        page.select_option('#cof-preset','tpb-dmtp-cof-pore-v1')
        wait("!COFdrawCrystalDesign.active && COFdrawExtent.visible && !COFdrawStudio.state.recipe?.assembly")
        save('returned-binary',2)
        assert not errors,errors
        print('Reaction SVG/PNG, 1/2/3 precursors, three pores, COF-300 depth, canvas switching and compact inspector passed.')
        browser.close()


if __name__=='__main__':
    main()
