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
            button=page.locator('#exports').get_by_role('button',name=ext.upper(),exact=True)
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
        page.locator('#measure-pore').click()
        wait("Boolean(COFdrawMeasure.active) && COFdrawScheme.composed")
        save('precursors-with-diameter',2)
        assert ET.parse(args.output/'precursors-with-diameter.svg').find(".//*[@data-cofdraw-pore-measurement='1']") is not None
        page.click('#figure-name')
        page.click('#figure-name-formula')
        assert page.locator('#figure-name-text').input_value()
        page.fill('#figure-name-text','TFB + TAB');page.fill('#figure-name-subscript','imine')
        page.click('#figure-name-apply')
        save('named-and-measured',2)
        doc=ET.parse(args.output/'named-and-measured.svg')
        name=doc.find(".//*[@data-cof-name='1']")
        assert name is not None and 'TFB + TAB' in ''.join(name.itertext())
        assert doc.find(".//*[@data-cofdraw-pore-measurement='1']") is not None
        page.get_by_role('button',name='Edit the name printed on the figure').dblclick()
        assert page.locator('#figure-name-text').input_value()=='TFB + TAB'
        page.fill('#figure-name-text','Edited framework')
        page.locator('#figure-name-text').press('Enter')
        with page.expect_download() as project:
            page.click('#save-project')
        filename=args.output/'named-project.json';project.value.save_as(filename)
        page.set_input_files('#project-file',str(filename))
        wait("COFdrawScheme.composed && !COFdrawStudio.state.busy && !COFdrawStudio.state.dirty")
        save('restored-name',2)
        assert 'Edited framework' in ''.join(ET.parse(args.output/'restored-name.svg').find(".//*[@data-cof-name='1']").itertext())
        page.click('#figure-name');page.click('#figure-name-remove')

        page.get_by_role('combobox',name='Drawing extent',exact=True).select_option('three-pores')
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
        wait("window.COFdrawCrystalDesign?.report && document.querySelector('#framework-image')?.naturalWidth>0")
        save('diamond-scheme',2,crystal=True)
        page.locator('#appearance-options summary').first.click()
        page.locator('#style-list [data-style="reticular_blocks"]').click()
        page.locator('#appearance-options summary').first.click()
        save('diamond-depth',2,crystal=True)
        page.locator('#open-science').click()
        wait("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        page.locator('#close-science').click()
        assert page.locator('#canvas-sheet.showing-scheme .reaction-figure').is_visible()
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
