from weasyprint import HTML, CSS

css = CSS(string=f'''
          @page {{
          size: {snakemake.params['width']}mm {snakemake.params['height']}mm;
          margin: 1in;
          }}

          body {{
          width: 100% !important;
          max-width: none !important;
          padding: 13px;
          }}

          .wide_table {{
          width: 100%; 
          max-width: 100%; 
          table-layout: fixed; 
          word-wrap: break-word;
          overflow-wrap: break-word;
          }}

          .wide_table th, .wide_table td {{
          white-space: normal;
          text-align: center;
          }}

          .cwm img {{
          max-width: 100px;
          height: auto;
          }}

          .distplot img {{
          max-width: 400px !important;
          height: auto;
          }}
''')

HTML(snakemake.input['html']).write_pdf(snakemake.output['pdf'], stylesheets=[css])
