"""
Configuration module for API keys and settings
Migrated from file-based to environment variable-based configuration
"""

import os
from pathlib import Path

# Try to load .env file if python-dotenv is available
try:
    from dotenv import load_dotenv
    # Load .env file from project root
    env_path = Path(__file__).parent.parent / ".env"
    if env_path.exists():
        load_dotenv(env_path)
except ImportError:
    pass  # python-dotenv not installed, use system environment variables


def get_deepseek_api_key() -> str:
    """
    Get DeepSeek API key from environment variable
    
    Returns:
        str: DeepSeek API key
        
    Raises:
        ValueError: If API key is not found
    """
    api_key = os.getenv("DEEPSEEK_API_KEY")
    if not api_key:
        raise ValueError(
            "DEEPSEEK_API_KEY not found in environment variables.\n"
            "Please set it using: export DEEPSEEK_API_KEY='your-key-here'"
        )
    return api_key


def get_grok_api_key() -> str:
    """
    Get Grok (xAI) API key from environment variable
    
    Returns:
        str: Grok API key
        
    Raises:
        ValueError: If API key is not found
    """
    api_key = os.getenv("GROK_API_KEY")
    if not api_key:
        raise ValueError(
            "GROK_API_KEY not found in environment variables.\n"
            "Please set it using: export GROK_API_KEY='your-key-here'\n"
            "Get your key at: https://console.x.ai/"
        )
    return api_key


def get_tavily_api_key() -> str:
    """
    Get Tavily API key from environment variable
    
    Returns:
        str: Tavily API key
        
    Raises:
        ValueError: If API key is not found
    """
    api_key = os.getenv("TAVILY_API_KEY")
    if not api_key:
        raise ValueError(
            "TAVILY_API_KEY not found in environment variables.\n"
            "Please set it using: export TAVILY_API_KEY='your-key-here'"
        )
    return api_key


def get_pubmed_email() -> str:
    """
    Get PubMed email from environment variable or settings file
    
    Returns:
        str: Email for PubMed API
    """
    # Try environment variable first
    email = os.getenv("PUBMED_EMAIL")
    if email:
        return email
    
    # Fallback to settings file
    try:
        import json
        settings_file = Path(__file__).parent.parent / "setting.json"
        with open(settings_file, "r") as f:
            settings = json.load(f)
        return settings.get("email", "your.email@example.com")
    except:
        return "your.email@example.com"


# For backward compatibility - try environment variable first, then file
def get_api_key_legacy() -> str:
    """
    Legacy function for backward compatibility
    Tries environment variable first, then falls back to file
    
    Returns:
        str: DeepSeek API key
    """
    # Try environment variable
    api_key = os.getenv("DEEPSEEK_API_KEY")
    if api_key:
        return api_key
    
    # Fallback to file (deprecated)
    try:
        api_key_file = Path(__file__).parent.parent / "api_key"
        with open(api_key_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        raise ValueError(
            "API key not found. Please set DEEPSEEK_API_KEY environment variable:\n"
            "export DEEPSEEK_API_KEY='your-key-here'"
        )
