# PendulumSim

A Julia project with AI agent integration via MCPRepl.

## Quick Start

```bash
cd PendulumSim
./repl
```

The MCP server will start automatically when Julia launches!

## Project Structure

```
PendulumSim/
├── src/                # Source code
├── test/               # Test suite
├── .mcprepl/           # Security configuration (git-ignored)
├── .vscode/            # VS Code MCP config
└── AGENTS.md           # AI agent guidelines
```

## Security

**Mode**: `lax` | **Port**: `3076` | **Auth**: None (localhost only)

To change security mode:

```julia
using MCPRepl
MCPRepl.setup()
```


## AI Agent Integration

**For AI agents**: See [AGENTS.md](AGENTS.md) for detailed guidelines.

**VS Code**: Open project, start Julia REPL  
**Claude Desktop**: Config in `.mcp.json`  
**Gemini**: Configured in `~/.gemini/settings.json`

## Troubleshooting

**Port in use?** Override with: `JULIA_MCP_PORT=3001 julia --project=.`

**Auth fails?** Check API key: `cat .env`

**Server won't start?** Restart Julia or check port: `lsof -i :3076`

## License

See LICENSE file
