# Publishing the VCFX Documentation

This document explains how to publish and maintain the VCFX documentation website.

## Overview of the Documentation System

The VCFX documentation is built using [MkDocs](https://www.mkdocs.org/) with the [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme. The system includes:

1. **Markdown Documentation Files**: Located in the `docs/` directory
2. **Configuration**: The `mkdocs.yml` file defines the site structure and navigation
3. **GitHub Actions Workflow**: Automated build and deployment process in `.github/workflows/docs.yml`

## Testing Locally

Before pushing any changes, you should test the documentation site locally:

1. Install the required dependencies:
   ```bash
   pip install mkdocs-material pymdown-extensions
   ```

2. Start the local development server:
   ```bash
   mkdocs serve
   ```

3. View the site in your browser at http://localhost:8000

4. Make changes to the documentation files and see them updated in real-time

## Publishing Process

### Automated Publishing (Recommended)

The documentation is automatically built and published when changes are pushed to the `main` branch. The GitHub Actions workflow will:

1. Build the documentation site
2. Deploy it to GitHub Pages

You can monitor the deployment progress in the "Actions" tab of the GitHub repository.

### Manual Publishing

If you need to publish the documentation manually:

1. Build the site:
   ```bash
   mkdocs build
   ```

2. Deploy to GitHub Pages:
   ```bash
   mkdocs gh-deploy --force
   ```

## Adding New Documentation

To add documentation for a new tool:

1. Create a new markdown file in the `docs/` directory named `VCFX_tool_name.md`
2. Add the tool to the appropriate section in the `nav` section of `mkdocs.yml`
3. Update any relevant index or overview pages to include the new tool

## Updating Existing Documentation

When updating existing documentation:

1. Edit the relevant Markdown files in the `docs/` directory
2. Test locally using `mkdocs serve`
3. Commit and push your changes to the repository

## Documentation Style Guidelines

To maintain consistent documentation:

1. Follow the established structure for tool documentation:
   - Overview
   - Usage
   - Options
   - Description
   - Examples
   - Limitations

2. Use consistent formatting for:
   - Command examples (code blocks with bash syntax highlighting)
   - Option tables (consistent column structure)
   - Section hierarchy (heading levels)

3. Include practical examples that demonstrate common use cases

## Troubleshooting

### Common Issues

**Documentation Not Updating**:
- Check the GitHub Actions workflow for errors
- Ensure the file is included in the navigation in `mkdocs.yml`
- Verify that the file paths are correct

**Local Preview Problems**:
- Make sure MkDocs and all required extensions are installed
- Check for Markdown syntax errors
- Restart the `mkdocs serve` command

**Deployment Failures**:
- Check GitHub Actions logs for specific errors
- Ensure GitHub Pages is configured correctly in the repository settings

## Resources

- [MkDocs Documentation](https://www.mkdocs.org/)
- [Material for MkDocs Documentation](https://squidfunk.github.io/mkdocs-material/)
- [Markdown Guide](https://www.markdownguide.org/)
- [GitHub Pages Documentation](https://docs.github.com/en/pages) 